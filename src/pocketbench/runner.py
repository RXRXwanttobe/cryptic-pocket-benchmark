"""Pinned inference adapters. These call upstream tools rather than reimplementing them."""
import collections,csv,dataclasses,gzip,hashlib,importlib.metadata,json,os,re,subprocess,time
from pathlib import Path
from .io import read,write,sha256

def preflight(root,models,p2rank=None,model_dir=None,allow_cpu=False):
    provenance=read(root/"config/original-provenance.json")
    runtime={"versions":{},"source_hashes":{},"allow_cpu":allow_cpu}
    if "p2rank" in models:
        if p2rank is None or not Path(p2rank).is_file():
            raise ValueError("Pass --p2rank /path/to/p2rank_2.5.1/prank")
        model=Path(p2rank).resolve().parent/"models/default/model.zst"
        if not model.exists() or sha256(model)!=provenance["p2rank_model_sha256"]:
            raise ValueError("P2Rank default model does not match the archived 2.5.1 model")
        version=subprocess.run([str(Path(p2rank).resolve()),"--version"],capture_output=True,text=True,timeout=60)
        runtime["p2rank_version_output"]=(version.stdout+version.stderr).strip()
        if "2.5.1" not in runtime["p2rank_version_output"]:
            raise ValueError("P2Rank 2.5.1 is required")
    if any(m!="p2rank" for m in models):
        import lacuna.pockets.clusterer as module
        if importlib.metadata.version("lacuna-pockets")!="1.2.0":
            raise ValueError("This protocol requires lacuna-pockets==1.2.0")
        for name,expected in provenance["source_hashes"].items():
            actual=sha256(Path(module.__file__).with_name(name))
            if actual!=expected:raise ValueError(f"Lacuna file differs from the archived execution: {name}")
            runtime["source_hashes"][name]=actual
        for name in ["lacuna-pockets","numpy","scipy","biopython"]:
            runtime["versions"][name]=importlib.metadata.version(name)
    if "learned_plm" in models:
        import torch
        import lacuna.pockets.plm as plm
        if not torch.cuda.is_available() and not allow_cpu:
            raise ValueError("Original PLM inference used CUDA FP16. Use --allow-cpu explicitly for a different numerical path")
        if model_dir is None:
            raise ValueError("Pass --model-dir after running download-model")
        model_dir=Path(model_dir).resolve()
        for name,info in provenance["esm"]["files"].items():
            path=model_dir/name
            if not path.exists() or sha256(path)!=info["sha256"]:
                raise ValueError(f"ESM checkpoint is missing or differs: {path}")
        plm.MODEL_NAME=str(model_dir)
        runtime.update(torch=torch.__version__,cuda=torch.version.cuda,
                       device=torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu",
                       precision="fp16" if torch.cuda.is_available() else "fp32",
                       esm_revision=provenance["esm"]["revision"])
        # The flag permits CPU when CUDA is absent; it does not disable an available GPU.
    return runtime

def lacuna_predictions(clusters):
    return [{"rank":c.rank,"residues":sorted({re.fullmatch(r"[A-Z]{3}(\d+):A",s)[1] for s in c.lining_residues},key=int)} for c in clusters]

def candidate_identity(c):
    data=c.to_dict()
    return json.dumps({k:v for k,v in data.items() if k not in ("rank","model_score","geometry_features","plm_features")},sort_keys=True)

def run_target(root,target,directory,models,p2rank,runtime):
    directory.mkdir(parents=True,exist_ok=True);input_file=root/target["input"]
    if sha256(input_file)!=target["input_sha256"]:raise ValueError(f"Input checksum mismatch: {target['id']}")
    config={"input_sha256":target["input_sha256"],"models":models,"runtime":runtime,"nma_conformers":20,"seed":20261023,"adapter_sha256":sha256(Path(__file__))}
    fingerprint=hashlib.sha256(json.dumps(config,sort_keys=True).encode()).hexdigest()
    status_file=directory/"predictions.json"
    if status_file.exists():
        existing=read(status_file)
        if existing.get("config_sha256")!=fingerprint:raise ValueError(f"Run configuration changed; use a new output folder: {directory}")
        if all(existing["status"][m]=="completed" for m in models):return existing
    record={"id":target["id"],"config_sha256":fingerprint,"models":{},"status":{},"errors":{},"timing_seconds":{}}
    if "p2rank" in models:
        started=time.monotonic()
        try:
            output=directory/"p2rank";output.mkdir(exist_ok=True)
            with (directory/"p2rank.log").open("w") as log:
                subprocess.run([str(Path(p2rank).resolve()),"predict","-f",str(input_file),"-o",str(output),"-threads","4","-visualizations","0"],stdout=log,stderr=subprocess.STDOUT,check=True,timeout=1800)
            files=list(output.glob("*_predictions.csv"))
            if len(files)!=1:raise ValueError("Expected one P2Rank predictions CSV")
            pockets=[]
            with files[0].open() as f:
                for row in csv.DictReader(f):
                    row={k.strip():v.strip() for k,v in row.items()}
                    pockets.append({"rank":int(row["rank"]),"residues":sorted({re.fullmatch(r"A_(\d+)",s)[1] for s in row["residue_ids"].split()},key=int)})
            record["models"]["p2rank"]=pockets;record["status"]["p2rank"]="completed"
        except Exception as exc:
            record["models"]["p2rank"]=[];record["status"]["p2rank"]="failed";record["errors"]["p2rank"]=str(exc)
        record["timing_seconds"]["p2rank"]=time.monotonic()-started
    lacuna_models=[m for m in models if m!="p2rank"]
    if lacuna_models:
        started=time.monotonic()
        try:
            from lacuna.io.structure import load_structure,coords_array
            from lacuna.ensemble.nma_backend import NMABackend
            from lacuna.pockets.detector import detect_pockets
            from lacuna.pockets.clusterer import cluster_pockets
            structure=load_structure(input_file)
            if [r.seq_num for r in structure.residues]!=list(range(1,len(structure.residues)+1)):
                raise ValueError("Use the bundled canonically numbered input")
            conformers=[coords_array(structure)]+list(NMABackend(seed=20261023).generate(input_file,20))
            candidates=[detect_pockets(coords,structure) for coords in conformers]
            def convert(x):
                if hasattr(x,"tolist"):return x.tolist()
                if hasattr(x,"item"):return x.item()
                raise TypeError(type(x))
            with gzip.open(directory/"candidate-pools.json.gz","wt") as f:
                json.dump([[dataclasses.asdict(p) for p in pockets] for pockets in candidates],f,default=convert)
            basic=cluster_pockets(candidates,21,rank_by="learned")
            record["timing_seconds"]["shared_geometry"]=time.monotonic()-started
            if "learned" in models:
                record["models"]["learned"]=lacuna_predictions(basic);record["status"]["learned"]="completed"
            if "learned_plm" in models:
                started=time.monotonic()
                try:
                    from lacuna.pockets import plm
                    probabilities=plm.residue_probabilities(structure)
                    if len(probabilities)!=len(structure.residues):raise ValueError("ESM-to-residue mapping incomplete")
                    augmented=cluster_pockets(candidates,21,rank_by="learned-plm",plm_residue_probs=probabilities)
                    if collections.Counter(map(candidate_identity,basic))!=collections.Counter(map(candidate_identity,augmented)):
                        raise ValueError("The complete candidate pool changed during PLM reranking")
                    write(directory/"plm-residue-scores.json",probabilities)
                    record["models"]["learned_plm"]=lacuna_predictions(augmented);record["status"]["learned_plm"]="completed"
                    record["shared_lacuna_candidates_verified"]=True
                except Exception as exc:
                    record["models"]["learned_plm"]=[];record["status"]["learned_plm"]="failed";record["errors"]["learned_plm"]=str(exc)
                record["timing_seconds"]["plm"]=time.monotonic()-started
        except Exception as exc:
            for m in lacuna_models:
                record["models"].setdefault(m,[]);record["status"].setdefault(m,"failed");record["errors"].setdefault(m,str(exc))
    write(status_file,record)
    return record
