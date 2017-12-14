import random
import string
from storage import DictListStorage, DictSetStorage
_integration_precision = 0.001

def _random_name(length):
    # For use with Redis, we return bytes
    return ''.join(random.choice(string.ascii_lowercase)
                   for _ in range(length)).encode('utf8')

def _integration(f, a, b):
    p = _integration_precision
    area = 0.0
    x = a
    while x < b:
        area += f(x + 0.5 * p) * p
        x += p
    return area, None


try:
    from scipy.integrate import quad as integrate
except ImportError:
    # For when no scipy installed
    integrate = _integration


def _false_positive_probability(threshold, b, r):
    _probability = lambda s: 1 - (1 - s ** float(r)) ** float(b)
    a, err = integrate(_probability, 0.0, threshold)
    return a


def _false_negative_probability(threshold, b, r):
    _probability = lambda s: 1 - (1 - (1 - s ** float(r)) ** float(b))
    a, err = integrate(_probability, threshold, 1.0)
    return a


def _optimal_param(threshold, num_perm, false_positive_weight,
                   false_negative_weight):
    '''
    Compute the optimal `MinHashLSH` parameter that minimizes the weighted sum
    of probabilities of false positive and false negative.
    '''
    min_error = float("inf")
    opt = (0, 0)
    best_fp = 1
    best_fn = 1
    for b in range(1, num_perm + 1):
        max_r = int(num_perm / b)
        for r in range(1, max_r + 1):
            fp = _false_positive_probability(threshold, b, r)
            fn = _false_negative_probability(threshold, b, r)
            best_fn = min(fn, best_fn)
            best_fp = min(fp, best_fp)
            error = fp * false_positive_weight + fn * false_negative_weight
            if error < min_error:
                min_error = error
                opt = (b, r)
    print("optimal number of band is {}, and each band have {} hashing".format(opt[0], opt[1]))
    print("optimal FN: {} and optimal FP: {}".format(best_fn, best_fp))
    return opt


class BioLSH(object):
    """
    Class for Locality Sensitive Hashing for DNA sequences
    """

    def __init__(self, threshold=0.9, num_perm=128):
        """
        create BioLSH object
        :param threshold: The Jaccard similarity threshold between 0.0 and
            1.0. The initialized MinHash LSH will be optimized for the threshold by
            minizing the false positive and false negative.
        :param num_perm: The number of permutation functions used
            by the MinHash to be indexed.
        """
        if threshold > 1.0 or threshold < 0.0:
            raise ValueError("threshold must be in [0.0, 1.0]")
        if num_perm < 2:
            raise ValueError("Too few permutation functions")

        self.h = num_perm
        self.b, self.r = _optimal_param(threshold, num_perm,
                                        0.5, 0.5)
        self.hashtables = [ DictSetStorage() for i in range(self.b)]
        self.hashranges = [(i * self.r, (i + 1) * self.r) for i in range(self.b)]
        self.keys = DictListStorage()

    def insert(self, key, minhash):
        '''
        Insert a unique key to the index, together
        with a MinHash (or weighted MinHash) of the set referenced by
        the key.
        Args:
            key (hashable): The unique identifier of the set.
            minhash (datasketch.MinHash): The MinHash of the set.
        '''
        self._insert(key, minhash, check_duplication=True, buffer=False)

    def _insert(self, key, minhash, check_duplication=True, buffer=False):
        if len(minhash) != self.h:
            raise ValueError("Expecting minhash with length %d, got %d"
                             % (self.h, len(minhash)))
        if check_duplication and key in self.keys:
            raise ValueError("The given key already exists")
        Hs = [self._H(minhash.hashvalues[start:end])
              for start, end in self.hashranges]

        # print("{}'s LSH bucket is {}".format(key, Hs))

        self.keys.insert(key, *Hs, buffer=buffer)
        for H, hashtable in zip(Hs, self.hashtables):
            hashtable.insert(H, key, buffer=buffer)

    def query(self, minhash):
        '''
        Giving the MinHash of the query set, retrieve
        the keys that references sets with Jaccard
        similarities greater than the threshold.

        Args:
            minhash (datasketch.MinHash): The MinHash of the query set.
        Returns:
            `list` of keys.
        '''
        if len(minhash) != self.h:
            raise ValueError("Expecting minhash with length %d, got %d"
                             % (self.h, len(minhash)))
        candidates = set()
        for (start, end), hashtable in zip(self.hashranges, self.hashtables):
            H = self._H(minhash.hashvalues[start:end])
            print(H)
            for key in hashtable.get(H):
                candidates.add(key)
        else:
            return list(candidates)

    def cluster(self):
        """
        clustering the input based on if they on same busket after lsh
        :return: clusters: list of cluster
        """
        clusters = []
        for (start, end), hashtable in zip(self.hashranges, self.hashtables):
            for key in hashtable.keys():
                clusters.append(hashtable[key])

        return clusters

    def insertion_session(self):
        '''
        Create a context manager for fast insertion into this index.
        Returns:
            datasketch.lsh.MinHashLSHInsertionSession
        '''
        return MinHashLSHInsertionSession(self)

    def __contains__(self, key):
        '''
        Args:
            key (hashable): The unique identifier of a set.
        Returns:
            bool: True only if the key exists in the index.
        '''
        return key in self.keys

    def remove(self, key):
        '''
        Remove the key from the index.
        Args:
            key (hashable): The unique identifier of a set.
        '''
        if key not in self.keys:
            raise ValueError("The given key does not exist")
        for H, hashtable in zip(self.keys[key], self.hashtables):
            hashtable.remove_val(H, key)
            if not hashtable.get(H):
                hashtable.remove(H)
        self.keys.remove(key)

    def is_empty(self):
        '''
        Returns:
            bool: Check if the index is empty.
        '''
        return any(t.size() == 0 for t in self.hashtables)

    @staticmethod
    def _H(hs):
        """
        Not hashing the band but instaed swap all values to generate new sequences
        :param hs:
        :return:
        """
        return bytes(hs.byteswap().data)

    def _query_b(self, minhash, b):
        if len(minhash) != self.h:
            raise ValueError("Expecting minhash with length %d, got %d"
                             % (self.h, len(minhash)))
        if b > len(self.hashtables):
            raise ValueError("b must be less or equal to the number of hash tables")
        candidates = set()
        for (start, end), hashtable in zip(self.hashranges[:b], self.hashtables[:b]):
            H = self._H(minhash.hashvalues[start:end])
            if H in hashtable:
                for key in hashtable[H]:
                    candidates.add(key)
        else:
            return candidates

    def get_counts(self):
        '''
        Returns a list of length ``self.b`` with elements representing the
        number of keys stored under each bucket for the given permutation.
        '''
        counts = [
            hashtable.itemcounts() for hashtable in self.hashtables]
        return counts

    def get_subset_counts(self, *keys):
        '''
        Returns the bucket allocation counts (see :ref:`get_counts` above)
        restricted to the list of keys given.
        Args:
            keys (hashable) : the keys for which to get the bucket allocation
                counts
        '''

        key_set = list(set(keys))
        hashtables = [DictSetStorage() for _ in
                      range(self.b)]
        Hss = self.keys.getmany(*key_set)
        for key, Hs in zip(key_set, Hss):
            for H, hashtable in zip(Hs, hashtables):
                hashtable.insert(H, key)
        return [hashtable.itemcounts() for hashtable in hashtables]

class MinHashLSHInsertionSession:
    '''Context manager for batch insertion of documents into a MinHashLSH.
    '''

    def __init__(self, lsh):
        self.lsh = lsh

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.lsh.keys.empty_buffer()
        for hashtable in self.lsh.hashtables:
            hashtable.empty_buffer()

    def insert(self, key, minhash, check_duplication=True):
        '''
        Insert a unique key to the index, together
        with a MinHash (or weighted MinHash) of the set referenced by
        the key.
        Args:
            key (hashable): The unique identifier of the set.
            minhash (datasketch.MinHash): The MinHash of the set.
        '''
        self.lsh._insert(key, minhash, check_duplication=check_duplication,
                         buffer=True)


if __name__ == '__main__':
    from minhash import BioMinHash
    from DNA_hashing import DNAHash

    set1 = set(["ATTTTT", "AGGGCT", "TTAAAA"])
    m1 = BioMinHash(num_perm=4)
    for d in set1:
        m1.update(d.encode('utf8'))
    set2 = set(["ATTTTT", "AGGGCT", "TAAAAA"])
    m2 = BioMinHash(num_perm=4)
    for d in set2:
        m2.update(d.encode('utf8'))

    lsh = BioLSH(threshold=0.2, num_perm=4)
    lsh.insert("m1", m1)
    result = lsh.cluster()
    print(result)