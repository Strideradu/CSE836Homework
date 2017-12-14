import random, copy, struct
from hashlib import sha1
from DNA_hashing import DNAHash
import numpy as np
import scipy.spatial.distance

# The size of a hash value in number of bytes
hashvalue_byte_size = len(bytes(np.int64(42).data))

# http://en.wikipedia.org/wiki/Mersenne_prime
_mersenne_prime = (1 << 61) - 1
_max_hash = (1 << 32) - 1
_hash_range = (1 << 32)


def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh'))[::-1].encode('utf8')


class BioMinHash(object):
    """
    A class to generate minhash signature for given biology sequences
    """

    def __init__(self, num_perm=128, seed=1, hashobj="DNAHash"):
        """

        :param num_perm: Number of random permutation functions.
        :param seed: The random seed controls the set of random permutation functions generated for this MinHash.
        :param hashobj: The hash function used by this MinHash. It must implements the `digest()` method similar
            to hashlib_ hash functions, such as `hashlib.sha1`.
        """
        if num_perm > _hash_range:
            # Because 1) we don't want the size to be too large, and
            # 2) we are using 4 bytes to store the size value
            raise ValueError("Cannot have more than %d number of\
                            permutation functions" % _hash_range)

        self.seed = seed
        if hashobj == "DNAHash":
            self.hashobj = DNAHash
        elif hashobj == "sha1":
            self.hashobj = sha1
        # Initialize hash values
        self.hashvalues = self._init_hashvalues(num_perm)
        generator = np.random.RandomState(self.seed)
        # Create parameters for a random bijective permutation function
        # that maps a 32-bit hash value to another 32-bit hash value.
        # http://en.wikipedia.org/wiki/Universal_hashing
        self.permutations = np.array([(generator.randint(1, _mersenne_prime, dtype=np.uint64),
                                       generator.randint(0, _mersenne_prime, dtype=np.uint64))
                                      for _ in range(num_perm)], dtype=np.uint64).T
        if len(self) != len(self.permutations[0]):
            raise ValueError("Numbers of hash values and permutations mismatch")

    def _init_hashvalues(self, num_perm):
        return np.ones(num_perm, dtype=np.uint64) * _max_hash

    def _parse_hashvalues(self, hashvalues):
        return np.array(hashvalues, dtype=np.uint64)

    def update(self, b, rc=False):
        '''
        Update this MinHash with a new value.

        Args:
            b (bytes): The value of type `bytes`.

        Example:
            To update with a new string value:

            .. code-block:: python
                minhash.update("new value".encode('utf-8'))
        '''
        # print(self.hashobj(b).digest()[:4])
        if self.hashobj == DNAHash:
            hv = struct.unpack('<I', self.hashobj(b, rc=rc).digest()[:4])[0]
        else:
            hv = struct.unpack('<I', self.hashobj(revcomp(b)).digest()[:4])[0]
        # print(hv)
        a, b = self.permutations
        phv = np.bitwise_and((a * hv + b) % _mersenne_prime, np.uint64(_max_hash))
        self.hashvalues = np.minimum(phv, self.hashvalues)

    def __len__(self):
        '''
        Returns:
            int: The number of hash values.
        '''
        return len(self.hashvalues)

    def jaccard(self, other):
        '''Estimate the `Jaccard similarity`_ (resemblance) between the sets
        represented by this MinHash and the other.

        Args:
            other (datasketch.MinHash): The other MinHash.

        Returns:
            float: The Jaccard similarity, which is between 0.0 and 1.0.
        '''
        if other.seed != self.seed:
            raise ValueError("Cannot compute Jaccard given MinHash with\
                    different seeds")
        if len(self) != len(other):
            raise ValueError("Cannot compute Jaccard given MinHash with\
                    different numbers of permutation functions")
        # return np.float(np.count_nonzero(self.hashvalues == other.hashvalues)) / np.float(len(self))
        return len(np.intersect1d(self.hashvalues, other.hashvalues)) / len(
            np.union1d(self.hashvalues, other.hashvalues))

    def identity(self, other):
        if other.seed != self.seed:
            raise ValueError("Cannot compute identity percentage given MinHash with\
                    different seeds")
        if len(self) != len(other):
            raise ValueError("Cannot compute identity percentage given MinHash with\
                    different numbers of permutation functions")
        return np.float(np.count_nonzero(self.hashvalues == other.hashvalues)) / np.float(len(self))

    def hamming(self, other):
        if other.seed != self.seed:
            raise ValueError("Cannot compute Hamming distance given MinHash with\
                    different seeds")
        if len(self) != len(other):
            raise ValueError("Cannot compute Hamming distance given MinHash with\
                    different numbers of permutation functions")

        return scipy.spatial.distance.hamming(self.hashvalues, other.hashvalues)

    def output(self):
        print("Frist hashing function is a = {} and b = {}, divisor is {}".format(self.permutations[0][0],
                                                                                  self.permutations[1][0],
                                                                                  _mersenne_prime))
        print("Signature vector is {}".format(self.hashvalues))
        print()

    def set_hashvalues(self, hashvalues):
        # only for debug purpose
        self.hashvalues = hashvalues


if __name__ == '__main__':
    set1 = set(["ATTTTT", "AGGGCT", "TTAAAA"])
    m1 = BioMinHash(num_perm=128)
    for d in set1:
        m1.update(d.encode('utf8'))

    set2 = set(["ATTTTT", "AGGGCT", "TAAAAA"])
    m2 = BioMinHash(num_perm=128)
    for d in set2:
        m2.update(d.encode('utf8'))

    print(m2.jaccard(m1))
