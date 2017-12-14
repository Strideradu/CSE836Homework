DNA_dict = {65:0, 67:1, 71:2, 84:3, "A":0, "C": 1, "G":2, "T":3}
RC_DNA_dict = {65:3, 67:2, 71:1, 84:0, "A":3, "C": 2, "G":1, "T":0}

class DNAHash(object):
    """
    Simple hashing function for DNA
    """
    def __init__(self, kmer, rc = False):
        """
        Convert input DNA kmer to 2k bit integer
        :param kmer:
        """
        self.k = len(kmer)
        self.int_value = 0
        if rc:
            for i, char in enumerate(kmer):
                self.int_value+=DNA_dict[char]*4**(i)
        else:
            for i, char in enumerate(kmer):
                self.int_value+=RC_DNA_dict[char]*4**(self.k - i - 1)
        self.invert_value = self._invert_hash(p = 2*self.k)

    def _invert_hash(self, p):
        """
        An invert hashing method from https://naml.us/post/inverse-of-a-hash-function/
        :param p:
        :return:
        """
        m = 2**p -1
        x = (~self.int_value + (self.int_value<< 21))&m
        x = x ^ x>>24
        x = (x + (x<<3) + (x<<8)) & m
        x = x ^ x>>14
        x = (x + (x << 2) + (x << 4)) & m
        x = x ^ x>>28
        x = (x+(x<<31)) & m
        return x

    def digest(self):
        # return (self.invert_value).to_bytes((2 * self.k + 7) //8, byteorder='big')
        return (self.invert_value).to_bytes(4, byteorder='big')

if __name__ == '__main__':

    test = DNAHash("ATCGCCCGCCCCGGGG")
    print(test.int_value)
    print(test.invert_value)
    print(test.digest())
    print(int.from_bytes(test.digest(), byteorder='big'))

