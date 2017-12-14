import collections


class FmCheckpoints(object):
    ''' Manages rank checkpoints and handles rank queries, which are
        O(1) time, with the checkpoints taking O(m) space, where m is
        length of text. '''

    def __init__(self, bw, cpIval=4):
        ''' Scan BWT, creating periodic checkpoints as we go '''
        self.cps = {}  # checkpoints
        self.cpIval = cpIval  # spacing between checkpoints
        tally = {}  # tally so far
        # Create an entry in tally dictionary and checkpoint map for
        # each distinct character in text
        for c in bw:
            if c not in tally:
                tally[c] = 0
                self.cps[c] = []
        # Now build the checkpoints
        for i in range(0, len(bw)):
            tally[bw[i]] += 1  # up to *and including*
            if (i % cpIval) == 0:
                for c in tally.keys():
                    self.cps[c].append(tally[c])

    def rank(self, bw, c, row):
        ''' Return # c's there are in bw up to and including row '''
        if row < 0 or c not in self.cps:
            return 0
        i, nocc = row, 0
        # Always walk to left (up) when calculating rank
        while (i % self.cpIval) != 0:
            if bw[i] == c:
                nocc += 1
            i -= 1
        return self.cps[c][i // self.cpIval] + nocc


class FMindex(object):
    def __init__(self, texts, ids):
        self.ids = {}
        s = []
        text = ''
        for i, t in enumerate(texts):

            if t[-1] != '$':
                t += '$'
            idx = len(self.ids)
            self.ids[idx] = ids[i]
            text += t
            s += [idx] * len(t)
        self.text = text
        self.text_len = len(text)
        self.suffix_array, self.sa_id = self._sa(texts)
        self.bwt, _ = self._bwt_and_dollar(texts, self.suffix_array, self.sa_id)

        self.cps = FmCheckpoints(self.bwt, 16)
        self.c = self.calculate_c_table()

    def _sa(self, s):
        ''' Given T return suffix array SA(T).  Uses "sorted"
            function for simplicity, which is probably very slow. '''
        satups = []
        for i, t in enumerate(s):
            if t[-1] != '$':
                t += '$'
            for j in range(0, len(t)):
                satups.append((t[j:], j, i))
        satups = sorted(satups)
        return list(map(lambda x: x[1], satups)), list(map(lambda x: x[2], satups))  # extract, return just offsets

    def _bwt(self, text, sa):
        return [text[i - 1] for i in sa]

    def _bwt_and_dollar(self, texts, sa, sa_id):
        ''' Given T, returns BWT(T) by way of the suffix array. '''
        bw = []
        dollarRow = {}
        for i, si in enumerate(sa):
            if si == 0:
                bw.append('$')
            else:
                bw.append(texts[sa_id[i]][si - 1])

        return (bw, dollarRow)

    def calculate_c_table(self):
        """returns a dict containing the counts of occurences of
        lexicographically lower valued characters for a given character"""
        tots = dict()
        for c in self.bwt:
            tots[c] = tots.get(c, 0) + 1
        # Calculate concise representation of first column
        first = {}
        totc = 0
        for c, count in sorted(tots.items()):
            first[c] = totc
            totc += count

        return first

    def count(self, c):
        ''' Return number of occurrences of characters < c '''
        if c not in self.c:
            # (Unusual) case where c does not occur in text
            for cc in sorted(self.c.keys()):
                if c < cc: return self.c[cc]
            return self.c[cc]
        else:
            return self.c[c]

    def range(self, p):
        ''' Return range of BWM rows having p as a prefix '''
        l, r = 0, self.text_len - 1  # closed (inclusive) interval
        for i in range(len(p) - 1, -1, -1):  # from right to left
            l = self.cps.rank(self.bwt, p[i], l - 1) + self.count(p[i])
            r = self.cps.rank(self.bwt, p[i], r) + self.count(p[i]) - 1
            if r < l:
                break
        return l, r + 1

    def update_range(self, p_i, l, r):
        l = self.cps.rank(self.bwt, p_i, l - 1) + self.count(p_i)
        r = self.cps.rank(self.bwt, p_i, r) + self.count(p_i) - 1
        return l, r

    def find_offset(self, row):
        ''' Given BWM row, return its offset w/r/t T '''

        def stepLeft(row):
            ''' Step left according to character in given BWT row '''
            c = self.bwt[row]
            return self.cps.rank(self.bwt, c, row - 1) + self.count(c)

        nsteps = 0
        while row not in self.suffix_array:
            row = stepLeft(row)
            nsteps += 1
        return self.suffix_array[row] + nsteps

    def hasSuffix(self, p):
        ''' Return true if and only if p is suffix of indexed text '''
        l, r = self.range(p)
        off = self.find_offset(l)
        return off, len(p)

    def find_sep(self, start, end):
        indices = self.suffix_array[start:end + 1]
        result = []
        for idx, i in enumerate(indices):
            if self.text[i - 1] == "$":
                result.append((i, self.sa_id[start + idx]))

        return result

    # in this function we backtrace the given pattern through the index
    # checking for prefixes at every step (indicated by '$')
    def find_overlaps(self, p, p_id, min_overlap_len):
        """returns all strings with prefix overlap over given min_overlap_len"""

        rs = reversed(p)
        start = 0
        end = self.text_len - 1
        overlap_size = 0
        all_overlaps = []
        for c in rs:

            start, end = self.update_range(c, start, end)
            # print(c, start, end)
            overlap_size += 1
            if start > end:
                break
            if min_overlap_len <= overlap_size < len(p):
                offs = self.find_sep(start, end)
                # print(offs)

                for off, read_idx in offs:
                    all_overlaps.append([p_id, str(len(p) - overlap_size ), str(len(p) - 1), self.ids[read_idx], str(0), str(overlap_size - 1), str(overlap_size)])

        return all_overlaps


if __name__ == '__main__':
    test_1 = "GTCAAAAAGTAGTGTGATTGGATGGCCTACTGTAAGGGAAAGAATGAGACGAGCTGAGCCAGCAGCAGATAGGGTGGGAGCAGCATCTCGAGACCTGGAAAAACATGGAGCAATCACAAGTAGCAATACAGCAGCTACCAATGCTGCTTGTGCCTGGCTAGAAGCACAAGAGGAGGAGGAGGTGGGTTTTCCAGTCACACCTCAGGTACCTTTAAGACCAATGACTTACAAGGCAGCTGTAGATCTTAGC"
    test_2 = "AAAAGGTCTATCTGGCATGGGTACCAGCACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGG"
    test_3 = 'ACACAAAGGAATTGGAGGAAATGAACAAGTAGATAAATTAGTCAGTGCTGGAATCAGGAAAGTACTATTTTTAGATGGAATAGATAAGGCCCAAGATGAACATGAGAAATATCACAGTAATTGGAGAGCAATGGCTAGTGATTTTAACCTGCCACCTGTAGTAGCAAAAGAAATAGTAGCCAGCTGTGATAAATGTCAGCTAAAAGGAGAAGCCATGCATGGACAAGTAGACTGTAGTCCAGGAATATGG'
    test1 = 'CCCGGGGTTTAC'
    test2 = 'CCCCCGGGGTTTA'
    test_index = FMindex([test_1, test_3], ['test1', 'test3'])
    test3 = 'TCCCCGGG'
    print(test_index.find_overlaps(test_2, 'test2', 200))
    #query_index = FMindex(test_1, "test_1")
    #print(test_index.find_overlaps(test_1, "test_1", 100))
