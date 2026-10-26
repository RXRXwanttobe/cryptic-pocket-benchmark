"""Command-line entry point; paths are relative to the chosen repository root."""
import argparse,csv,json,os,sys,tarfile
from pathlib import Path

for name in ["OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"]:os.environ[name]="4"
from .io import read,write,sha256,cohort,prepare_references,download
from .metrics import MODELS,summarize,score_target

def verify(root):
    targets=cohort(root)
    if len(targets)!=100 or len({t["pdb"] for t in targets})!=100:
        raise ValueError("Expected 100 distinct frozen target structures")
    predictions=read(root/"data/published/predictions.json")["targets"]
    metrics=read(root/"data/published/target_metrics.json")["targets"]
    if [t["id"] for t in targets]!=[t["id"] for t in predictions] or [t["id"] for t in targets]!=[t["id"] for t in metrics]:
        raise ValueError("Target ordering or membership differs")
    for t,p in zip(targets,predictions):
        if sha256(root/t["input"])!=t["input_sha256"]:raise ValueError(f"Input mismatch: {t['id']}")
        present=set(t["original_to_normalized_residue"].values())
        if present!=set(t["label_to_normalized_residue"].values()):raise ValueError("Residue mapping incomplete")
        for model in MODELS:
            scores=p["models"][model];ranks=[s["rank"] for s in scores]
            if ranks!=list(range(1,len(ranks)+1)):raise ValueError("Unexpected candidate ranks")
            if any(not set(s["residues"])<=present for s in scores):raise ValueError("Unmapped prediction residue")
    stats=summarize(metrics);expected=read(root/"data/published/statistics.json")
    if stats["totals"]!=expected["totals"] or stats["thresholds"]!=expected["thresholds"]:
        raise ValueError("Summary does not reproduce")
    for a,b in zip(stats["comparisons"],expected["comparisons"]):
        for key in a:
            if isinstance(a[key],float):
                if abs(a[key]-b[key])>1e-12:raise ValueError(f"Paired statistic differs: {key}")
            elif a[key]!=b[key]:raise ValueError(f"Paired statistic differs: {key}")
    print("PASS: 100 input checksums, residue mappings, 300 prediction sets, Top-k totals and paired statistics")
    print("This offline check does not independently rescore IoU. Run prepare then evaluate --published for that check.")

def evaluate(root,run_dir=None):
    reference_file=root/"data/reference.generated.json"
    if not reference_file.exists():raise ValueError("Run 'pocketbench prepare' to retrieve and map the external reference annotations")
    references=read(reference_file);targets=cohort(root);failures=[]
    published=run_dir is None
    if published:
        predictions=read(root/"data/published/predictions.json")["targets"];models=list(MODELS)
        destination=root/"runs/published-rescore"
    else:
        destination=run_dir;manifest=read(run_dir/"run.json");models=manifest["models"]
        targets=[t for t in targets if t["id"] in manifest["target_ids"]];predictions=[]
        for t in targets:
            path=run_dir/t["id"]/"predictions.json"
            if not path.exists():
                failures.append({"id":t["id"],"status":"missing"});predictions.append({"id":t["id"],"models":{m:[] for m in models}})
            else:
                p=read(path);predictions.append(p)
                if any(p["status"][m]!="completed" for m in models):failures.append({"id":t["id"],"status":p["status"],"errors":p["errors"]})
    by_id={p["id"]:p for p in predictions};rows=[]
    for t in targets:
        if sha256(root/t["input"])!=t["input_sha256"]:raise ValueError("Input changed before evaluation")
        row={k:t[k] for k in ["id","pdb","group","name"]};present=set(t["original_to_normalized_residue"].values())
        row["metrics"]={m:score_target(by_id[t["id"]]["models"][m],references[t["id"]],present) for m in models};rows.append(row)
    stats=summarize(rows,models)
    if published:
        frozen={r["id"]:r for r in read(root/"data/published/target_metrics.json")["targets"]}
        for row in rows:
            for m in models:
                for key,value in row["metrics"][m].items():
                    if value!=frozen[row["id"]]["metrics"][m][key]:raise ValueError(f"Rescored result differs: {row['id']} / {m} / {key}")
    write(destination/"target-metrics.json",{"targets":rows});write(destination/"statistics.json",stats);write(destination/"failures.json",failures)
    print(json.dumps({"n":len(rows),"failed_or_missing_targets":len(failures),"top135":{m:{str(k):stats["totals"][m][str(k)]["hits"] for k in [1,3,5]} for m in models}},indent=2))
    if published:print("PASS: all 300 prediction sets rescored against the pinned external reference")

def run(root,args):
    from .runner import preflight,run_target
    models=args.models.split(",")
    if not models or len(models)!=len(set(models)) or any(m not in MODELS for m in models):raise ValueError("Choose unique model names: p2rank,learned,learned_plm")
    models=[m for m in MODELS if m in models]
    selected=cohort(root)
    if args.target:
        wanted=set(args.target);selected=[t for t in selected if t["id"] in wanted]
        if {t["id"] for t in selected}!=wanted:raise ValueError("Unknown target ID")
    runtime=preflight(root,models,args.p2rank,args.model_dir,args.allow_cpu)
    run_dir=Path(args.out).resolve();manifest={"target_ids":[t["id"] for t in selected],"models":models,"runtime":runtime}
    if (run_dir/"run.json").exists() and read(run_dir/"run.json")!=manifest:raise ValueError("Output folder belongs to a different run; choose a fresh folder")
    write(run_dir/"run.json",manifest)
    for i,t in enumerate(selected,1):
        record=run_target(root,t,run_dir/t["id"],models,args.p2rank,runtime)
        print(f"[{i}/{len(selected)}] {t['id']} {t['pdb']}: {record['status']}",flush=True)
    print("Run complete. Evaluate with: pocketbench evaluate --run-dir "+str(run_dir))

def download_model(root,destination):
    provenance=read(root/"config/original-provenance.json")["esm"]
    destination=Path(destination).resolve()
    for name,info in provenance["files"].items():
        print("Checking/downloading "+name,flush=True)
        download(f"https://huggingface.co/{provenance['model']}/resolve/{provenance['revision']}/{name}",destination/name,info["sha256"])
    print("Verified pinned ESM-2 checkpoint: "+str(destination))

def install_p2rank(destination):
    destination=Path(destination).resolve();destination.mkdir(parents=True,exist_ok=True)
    archive=download("https://github.com/rdk/p2rank/releases/download/2.5.1/p2rank_2.5.1.tar.gz",destination/"p2rank_2.5.1.tar.gz","d243f2d9036ac053fefb9407b5fe1c85f4fe077c519fd975ac585e995feab274")
    if (destination/"p2rank_2.5.1").exists():raise ValueError("P2Rank directory already exists; use the existing installation")
    with tarfile.open(archive) as tar:
        for member in tar.getmembers():
            if not (destination/member.name).resolve().is_relative_to(destination):raise ValueError("Unsafe archive path")
        tar.extractall(destination,filter="data")
    print("Installed "+str(destination/"p2rank_2.5.1/prank"))

def main():
    parser=argparse.ArgumentParser(description="Reproduce a frozen comparison of P2Rank and two Lacuna rankers")
    parser.add_argument("--root",type=Path,default=Path.cwd(),help="Repository root; default: current directory")
    sub=parser.add_subparsers(dest="command",required=True)
    sub.add_parser("verify",help="Verify bundled inputs, stored results, and statistics offline")
    sub.add_parser("prepare",help="Explicitly download the pinned external CryptoBench reference annotations")
    p=sub.add_parser("evaluate",help="Rescore published predictions or a new run")
    p.add_argument("--published",action="store_true",help="Use the bundled predictions (default)")
    p.add_argument("--run-dir",type=Path)
    p=sub.add_parser("run",help="Run upstream predictors on the frozen input structures")
    p.add_argument("--models",default=",".join(MODELS));p.add_argument("--target",action="append",help="Target ID; repeat to choose multiple targets")
    p.add_argument("--p2rank",type=Path);p.add_argument("--model-dir",type=Path);p.add_argument("--allow-cpu",action="store_true")
    p.add_argument("--out",required=True,help="New output directory; matching completed targets can be resumed")
    p=sub.add_parser("download-model",help="Download the pinned approximately 2.6 GB ESM-2 encoder; verify SHA-256")
    p.add_argument("--out",default="models/esm2_t33_650M_UR50D")
    p=sub.add_parser("install-p2rank",help="Download and safely extract upstream P2Rank 2.5.1")
    p.add_argument("--out",default="tools")
    args=parser.parse_args();root=args.root.resolve()
    try:
        if args.command=="verify":verify(root)
        elif args.command=="prepare":print(f"Mapped references for {len(prepare_references(root))} targets")
        elif args.command=="evaluate":
            if args.published and args.run_dir:raise ValueError("Choose --published or --run-dir")
            evaluate(root,args.run_dir.resolve() if args.run_dir else None)
        elif args.command=="run":run(root,args)
        elif args.command=="download-model":download_model(root,args.out)
        elif args.command=="install-p2rank":install_p2rank(args.out)
    except Exception as exc:
        parser.exit(1,f"Error: {exc}\n")
