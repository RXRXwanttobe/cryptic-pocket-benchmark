"""Small, explicit file and download helpers."""
import hashlib,json,urllib.request
from pathlib import Path

def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))

def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+".tmp")
    temporary.write_text(json.dumps(value,indent=2,ensure_ascii=True)+"\n",encoding="utf-8")
    temporary.replace(path)

def sha256(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream,"sha256").hexdigest()

def download(url,path,expected_sha):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists():
        if sha256(path)!=expected_sha:
            raise ValueError(f"Existing file has an unexpected SHA-256: {path}")
        return path
    temporary=path.with_suffix(path.suffix+".partial")
    with urllib.request.urlopen(url,timeout=120) as source,temporary.open("wb") as destination:
        while chunk:=source.read(8*1024*1024):destination.write(chunk)
    if sha256(temporary)!=expected_sha:
        raise ValueError(f"Download checksum mismatch: {url}")
    temporary.replace(path)
    return path

def cohort(root):
    return read(root/"data/cohort.json")["targets"]

def prepare_references(root):
    source=read(root/"config/reference-source.json")
    file=download(source["url"],root/".cache/cryptobench-labels.json",source["sha256"])
    labels=read(file);references={}
    for target in cohort(root):
        raw=labels[target["pdb"].lower()]
        if {s.rsplit("_",1)[0] for s in raw}!={target["chain"]}:
            raise ValueError(f"Reference chain mismatch: {target['id']}")
        reference=sorted({target["label_to_normalized_residue"][s] for s in raw},key=int)
        digest=hashlib.sha256(json.dumps(reference).encode()).hexdigest()
        if len(reference)!=target["reference_count"] or digest!=target["reference_sha256"]:
            raise ValueError(f"Reference mapping differs from the frozen run: {target['id']}")
        references[target["id"]]=reference
    write(root/"data/reference.generated.json",references)
    return references
