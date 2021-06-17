
import os

from collections import namedtuple

from whoosh.fields import SchemaClass, TEXT, KEYWORD, ID
from whoosh.qparser import MultifieldParser
from whoosh import index


class SearchBarSchema(SchemaClass):
    seqid = ID(stored=True)
    locus_tag = KEYWORD(stored=True)
    gene = KEYWORD(stored=True)
    product = TEXT(stored=True)


SearchResult = namedtuple("SearchResult", ["seqid", "locus_tag", "gene", "product"])


class ChlamdbIndex:

    def new_index(name):
        chlamdb_index = ChlamdbIndex()
        chlamdb_index.index = index.create_in(name, SearchBarSchema)
        chlamdb_index.writer = chlamdb_index.index.writer()
        return chlamdb_index

    def use_index(name):
        chlamdb_index = ChlamdbIndex()
        chlamdb_index.index = index.open_dir(name)
        return chlamdb_index

    def add(self, **kwargs):
        self.writer.add_document(**kwargs)

    def search(self, user_query, limit=10):
        parser = MultifieldParser(["locus_tag", "gene", "product"], self.index.schema)
        query = parser.parse(user_query)

        for result in self.index.searcher().search(query, limit=limit):
            gene = result.get("gene", None)
            locus_tag = result.get("locus_tag", None)
            product = result.get("product", None)
            seqid = result.get("seqid", None)
            yield SearchResult(seqid=int(seqid), locus_tag=locus_tag, gene=gene, product=product)

    def done_adding(self):
        self.writer.commit(optimize=True)
