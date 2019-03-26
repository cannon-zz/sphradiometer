function run_transcode() {
	transcode -i "${1}" -g 800x800 -f 10 -Z 600x600 -w 2048 --use_rgb -x imlist,null -y xvid -z -o "${2}"
	python ../bin/make_subtitle.py <"${1}" >"${2/.avi/.srt}"
}

#ls power_td*.png | sort >images.list
#run_transcode images.list movie_td.avi
ls power_fd*.png | sort >images.list
run_transcode images.list movie_fd.avi
