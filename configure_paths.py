"""Create a separately configured copy without modifying the release archive."""
import argparse
from pathlib import Path
import shutil

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, required=True)
    for option in ('project', 'scratch', 'phenotypes', 'remote', 'micromamba', 'compute-host'):
        parser.add_argument('--' + option, required=True)
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    target = args.destination.resolve()
    if target == source or source in target.parents:
        parser.error('Choose a destination outside the release directory.')
    if target.exists():
        parser.error('The destination must not already exist.')
    mapping = {'/work/project': args.project, '/work/scratch': args.scratch,
               '/work/phenotypes': args.phenotypes, '/work/remote': args.remote,
               '/opt/micromamba': args.micromamba, 'compute-host': args.compute_host}
    for value in mapping.values():
        if any(char in value for char in ('"', "'", '`', '$', '\n', '\r', '\\')):
            parser.error('Use POSIX paths and a host alias without quotes or shell metacharacters.')
    shutil.copytree(source / 'gwas_frost', target / 'gwas_frost')
    for path in (target / 'gwas_frost/scripts').iterdir():
        content = path.read_text(encoding='utf-8')
        for old, new in mapping.items():
            content = content.replace(old, new)
        path.write_text(content, encoding='utf-8', newline='\n')
    print(target)

if __name__ == '__main__':
    main()
