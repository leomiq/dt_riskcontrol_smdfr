"""Dump the structure of a Simulink .slx file without MATLAB.

Reads the zipped XML of an .slx and prints, for every subsystem:
  - each block (SID, type, name, key parameters),
  - each block's input ports with the source block/port that drives them,
  - the code of every MATLAB Function block (Stateflow chart scripts).

Usage:
    python tools/slx_dump.py models/SMDFR_hrm_MCDET_2022b_control_CMEX.slx > hrm_dump.txt

Used in F1b (docs/monolithic_formulation_plan.md §8.6) to trace each equation of the
HRM replica to its block, as F0 did for the ROM.
"""
import html
import re
import sys
import zipfile
import xml.etree.ElementTree as ET

SKIP_PARAMS = {"Position", "ZOrder", "Ports", "BackgroundColor", "ForegroundColor",
               "ShowName", "HideAutomaticName", "LibraryVersion", "SourceBlock",
               "SourceType", "ContentPreviewEnabled", "RequestExecContextInheritance",
               "SFBlockType", "Permissions", "ErrorFcn", "PermitHierarchicalResolution",
               "TreatAsAtomicUnit", "MinAlgLoopOccurrences", "PropExecContextOutsideSubsystem",
               "SystemSampleTime", "RTWSystemCode", "RTWFcnNameOpts", "RTWFileNameOpts",
               "FunctionWithSeparateData", "IsSubsystemVirtual", "Variant", "BlockMirror",
               "BlockRotation", "NamePlacement", "FontName", "FontSize", "DropShadow",
               "DataTypeOverride", "MinMaxOverflowLogging", "DataTypeOverrideAppliesTo",
               "SFunctionDeploymentMode", "VariantControl", "CompiledIsActive"}


def read(z, name):
    return ET.fromstring(z.read(name))


def main(path):
    z = zipfile.ZipFile(path)
    names = z.namelist()
    systems = {n: read(z, n) for n in names if n.startswith("simulink/systems/system_")}

    # SID -> (system file, block element)
    blocks = {}
    for sysname, root in systems.items():
        for b in root.iter("Block"):
            blocks[b.get("SID")] = (sysname, b)

    # Which block owns each system file (subsystem hierarchy)
    owner = {}
    for sid, (sysname, b) in blocks.items():
        for s in b.iter("System"):
            ref = s.get("Ref")
            if ref:
                owner["simulink/systems/" + ref + ".xml"] = sid

    def fullname(sid):
        sysname, b = blocks[sid]
        parts = [b.get("Name")]
        while sysname in owner:
            sid2 = owner[sysname]
            sysname, b2 = blocks[sid2]
            parts.append(b2.get("Name"))
        return "/".join(reversed(parts))

    # Stateflow chart code keyed by chart name (= block path)
    charts = {}
    for n in names:
        if n.startswith("simulink/stateflow/chart_"):
            t = z.read(n).decode("utf8")
            nm = re.search(r'<P Name="name">([^<]*)', t)
            sc = re.search(r'<P Name="script">(.*?)</P>', t, re.S)
            if nm:
                charts[html.unescape(nm.group(1))] = html.unescape(sc.group(1)) if sc else ""

    def parse_endpoint(s):
        # "12#out:1" -> ("12", "out", "1")
        m = re.match(r"(\d+)#(\w+):(\w+)", s)
        return m.groups() if m else (s, "?", "?")

    for sysname in sorted(systems, key=lambda s: int(re.search(r"_(\d+)", s).group(1))
                          if re.search(r"_(\d+)\.xml", s) else 0):
        root = systems[sysname]
        title = fullname(owner[sysname]) if sysname in owner else "<root>"
        print("=" * 100)
        print(f"SYSTEM {title}   ({sysname})")
        # connections: dst -> src
        conns = {}

        def walk_line(el, src):
            for p in el.findall("P"):
                if p.get("Name") == "Src":
                    src = p.text
            for p in el.findall("P"):
                if p.get("Name") == "Dst":
                    conns.setdefault(p.text, []).append(src)
            for br in el.findall("Branch"):
                walk_line(br, src)

        for ln in root.findall("Line"):
            walk_line(ln, None)

        for b in root.findall("Block"):
            sid, bt, nm = b.get("SID"), b.get("BlockType"), b.get("Name")
            params = {p.get("Name"): (p.text or "") for p in b.findall("P")
                      if p.get("Name") not in SKIP_PARAMS}
            pstr = ", ".join(f"{k}={v}" for k, v in params.items() if len(v) < 120)
            print(f"\n[{sid}] {bt}: {nm!r}  {pstr}")
            ins = sorted((k for k in conns if k.startswith(sid + "#")),
                         key=lambda s: (s.split("#")[1].split(":")[0],
                                        int(s.split(":")[1]) if s.split(":")[1].isdigit() else 0))
            for k in ins:
                for src in conns[k]:
                    if src is None:
                        print(f"    {k.split('#')[1]:10s} <- ?")
                        continue
                    s_sid, _, s_port = parse_endpoint(src)
                    s_name = blocks[s_sid][1].get("Name") if s_sid in blocks else s_sid
                    print(f"    {k.split('#')[1]:10s} <- [{s_sid}] {s_name!r} out:{s_port}")

    print("=" * 100)
    print("MATLAB FUNCTION CODE")
    for nm in sorted(charts):
        print("-" * 100)
        print(nm)
        print(charts[nm])


if __name__ == "__main__":
    main(sys.argv[1])
