import requests


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


def download_txt_url(url):
        with requests.get(url, stream=True) as r:
            r.raise_for_status()
            for chunk in r.iter_content(chunk_size=8192, decode_unicode=True):
                # If you have chunk encoded response uncomment if
                # and set chunk_size parameter to None.
                # if chunk:
                yield chunk
