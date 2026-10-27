"""Stream a matching transcriptome and genome into a Salmon full-decoy reference.

No reference download or sequence renaming. Input FASTA IDs must be unique and
disjoint so Salmon cannot accidentally interpret a transcript as a decoy.
"""
from pathlib import Path
import gzip
import os


def prepare(transcripts, genome, directory):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    output = directory/'gentrome.fa'
    decoys = directory/'decoys.txt'
    seen = set()
    try:
        with output.with_suffix('.tmp').open('wb') as out, decoys.with_suffix('.tmp').open('wb') as names:
            for source, is_decoy in ((transcripts, False), (genome, True)):
                source = Path(source)
                with source.open('rb') as probe:
                    compressed = probe.read(2) == b'\x1f\x8b'
                count = 0
                with (gzip.open(source, 'rb') if compressed else source.open('rb')) as incoming:
                    for line in incoming:
                        if line.startswith(b'>'):
                            fields = line[1:].split()
                            if not fields or fields[0] in seen:
                                raise ValueError('FASTA 序列 ID 为空、重复，或基因组与转录本 ID 冲突')
                            identifier = fields[0]
                            seen.add(identifier)
                            count += 1
                            if is_decoy:
                                names.write(identifier+b'\n')
                        elif line.strip() and count == 0:
                            raise ValueError('参考文件不是有效 FASTA（缺少 > 序列头）')
                        out.write(line)
                    out.write(b'\n')
                if not count:
                    raise ValueError('参考 FASTA 没有序列')
        os.replace(output.with_suffix('.tmp'), output)
        os.replace(decoys.with_suffix('.tmp'), decoys)
    except BaseException:
        output.with_suffix('.tmp').unlink(missing_ok=True)
        decoys.with_suffix('.tmp').unlink(missing_ok=True)
        raise
    return output, decoys


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('transcripts')
    parser.add_argument('genome')
    parser.add_argument('directory')
    args = parser.parse_args()
    prepare(args.transcripts, args.genome, args.directory)
    print('已生成 full-decoy reference；使用匹配版本的原始 FASTA ID。')
