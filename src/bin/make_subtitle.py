import math
import sys

fps = 10.0


def frame_time(n):
	"""
	Time at which frame n is displayed, in seconds.
	"""
	return n / fps


def srt_time_string(t):
	"""
	Convert a time in seconds to a .srt style time string.
	"""
	h = int(t / 3600)
	t %= 3600
	m = int(t / 60)
	t %= 60
	s = int(t / 1)
	t %= 1
	t = int(1000 * t)
	return "%02d:%02d:%02d,%03d" % (h, m, s, t)


def subtitle(gmst):
	"""
	The subtitle string to display
	"""
	return "%.4f days" % (gmst / (2 * math.pi))


for n, filename in enumerate(sys.stdin):
	gmst = float(filename.strip()[9:-4])

	print n
	print srt_time_string(frame_time(n)) + " --> " + srt_time_string(frame_time(n + 1))
	print subtitle(gmst)
	print
