import unittest

class TopicTest(unittest.TestCase):
	def test_import(self):
		import topicpy
		self.assertEqual("1.1.6",topicpy.__version)
