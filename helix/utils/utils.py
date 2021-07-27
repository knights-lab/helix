import urllib


def read_fasta(fh):
    """
    :return: tuples of (title, seq)
    """
    title = None
    data = None
    for line in fh:
        if line[0] == ">":
            if title:
                yield title.strip(), data
            title = line[1:]
            data = ''
        else:
            data += line.strip()
    if not title:
        yield None
    yield title.strip(), data


def download_txt_url(path_to_file, url):
    with urllib.request.urlopen(url) as stream:
        CHUNK = 2 ** 14
        with open(path_to_file, 'wb') as outfile:
            while True:
                chunk = stream.read(CHUNK)
                if not chunk:
                    break
                outfile.write(chunk)
