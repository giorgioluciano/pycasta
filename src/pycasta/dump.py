from pathlib import Path
import os

ROOT = Path(__file__).resolve().parent
OUTFILE = ROOT / "project_dump.txt"

EXCLUDE_DIRS = {
    ".git", ".venv", "venv", "__pycache__",
    ".pytest_cache", ".mypy_cache",
    "build", "dist", "results", "result",
    "output", "data", ".idea", ".vscode",
}

EXCLUDE_FILES = {
    ".DS_Store",
}

TEXT_EXTS = {
    ".py", ".txt", ".md", ".toml", ".cfg", ".ini", ".yaml", ".yml",
    ".json", ".js", ".ts", ".sh", ".ps1", ".bat", ".cmd", ".rst",
}

def is_text_file(path: Path) -> bool:
    return path.suffix.lower() in TEXT_EXTS

def should_skip_dir(name: str) -> bool:
    return name in EXCLUDE_DIRS

def should_skip_file(name: str) -> bool:
    return name in EXCLUDE_FILES or name.endswith((".pyc", ".pyo", ".pyd"))

def build_tree(root: Path) -> str:
    lines = [root.name]

    def walk(dir_path: Path, prefix: str = ""):
        items = []
        for p in sorted(dir_path.iterdir(), key=lambda x: (not x.is_dir(), x.name.lower())):
            if p.is_dir() and should_skip_dir(p.name):
                continue
            if p.is_file() and should_skip_file(p.name):
                continue
            items.append(p)

        for i, p in enumerate(items):
            last = i == len(items) - 1
            branch = "└── " if last else "├── "
            lines.append(prefix + branch + p.name)
            if p.is_dir():
                walk(p, prefix + ("    " if last else "│   "))

    walk(root)
    return "\n".join(lines)

def collect_file_contents(root: Path) -> str:
    parts = []
    for dirpath, dirs, files in os.walk(root):
        dirs[:] = [d for d in dirs if not should_skip_dir(d)]
        for filename in sorted(files):
            if should_skip_file(filename):
                continue
            path = Path(dirpath) / filename
            if not is_text_file(path):
                continue
            try:
                content = path.read_text(encoding="utf-8", errors="replace")
            except Exception as e:
                content = f"[ERROR READING FILE: {e}]"
            parts.append(f"### FILE: {path.relative_to(root)}\n{content}\n")
    return "\n---\n\n".join(parts)

def main():
    tree = build_tree(ROOT)
    contents = collect_file_contents(ROOT)
    OUTFILE.write_text(
        f"ROOT: {ROOT}\n\nPROJECT TREE\n{tree}\n\n{'='*80}\n\nFILE CONTENTS\n{'='*80}\n\n{contents}\n",
        encoding="utf-8"
    )
    print(f"Written: {OUTFILE}")

if __name__ == "__main__":
    main()