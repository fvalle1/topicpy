import unittest

class TopicTest(unittest.TestCase):
	def test_import(self):
		import topicpy
		self.assertEqual("1.1.5",topicpy.__version)
