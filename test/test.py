import unittest
import topicpy

class TopicTest(unittest.TestCase):
	def test_load(self):
		self.assertEqual("1",topicpy.__version__[0])
		self.assertEqual("1",topicpy.__version__.split(".")[1])

	def test_hypergeom(self):
		import pandas as pd
		from topicpy.hypergeom import parameters_for_hypergeometric
		l1 = pd.Series(index=["ENSG00000000123", "ENSG00000000456", "ENSG00000000789", "ENSG00000000XXX"], data=["c1", "c1", "c1", "c2"], dtype=str)
		l2 = pd.Series(index=["ENSG00000000123", "ENSG00000000456", "ENSG00000000789"], data=["c1", "c1", "c1"], dtype=str)
		x, M, k, N, _ = parameters_for_hypergeometric(l1, l2)
		self.assertEqual(M, 3)
		self.assertEqual(k["c1"], 3)
		self.assertEqual(N["c1"],3)
		self.assertEqual(N["c2"],1)
		self.assertEqual(x.loc["c1","c1"], 3)

	def test_converter(self):
		from topicpy.converter import convert_list_to_sybmols, convert_list_to_ensg
		name = convert_list_to_sybmols(["ENSG00000159763"])
		self.assertEqual(name[0], "PIP")
		ensg = convert_list_to_ensg(["PIP"])
		self.assertEqual(ensg[0], "ENSG00000159763")

	def test_go(self):
		import pandas as pd
		from topicpy.geneontology import get_ontology_df
		df = get_ontology_df(["ENSG00000159763"])
		self.assertEqual(type(df), pd.DataFrame)

	def test_hsbmpy(self):
		from topicpy.hsbmpy import clusteranalysis
		clusteranalysis("test/fake",["SMTS","SMTSD"])
		self.assertEqual(True, True)

if __name__ == '__main__':
    unittest.main()
