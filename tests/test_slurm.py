import unittest


class TestSlurmMethods(unittest.TestCase):

	def test_loadDefault(self):
		opts = slurm.loadDefault( pkg_resources.resource_filename('ScanpyAutoAnalyzer',join( 'data', 'testScript.sh' )) )
		exp = {
			'n' : "10",
			'N' : "1",
			't' : "12:00:00",
			'A' : "lsens2018-3-3",
			'p' : "dell",
			'J' : "NA"
		}
		self.assertEqual( opts, exp )




if __name__ == '__main__':
    unittest.main()