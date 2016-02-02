
all:
	p4 edit *
	python ../test_bam_validation.py --make-files
	p4 revert -a ./...

*:
	make all
