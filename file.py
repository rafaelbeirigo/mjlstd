import pickle


def load(filename):
    """Loads persisted data."""
    try:
        with open(filename, 'rb') as f:
            return pickle.load(f)
    except Exception:
        return None


def save(data, filename):
    """Persists data."""
    with open(filename, 'wb') as f:
        pickle.dump(data, f, pickle.HIGHEST_PROTOCOL)


def saverep(data, name, r):
    save(data, 'pkl/' + name + '-{:03d}.pkl'.format(r + 1))


def loadrep(name, r):
    return load('pkl/' + name + '-{:03d}.pkl'.format(r + 1))
