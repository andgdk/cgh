run:
	scilab -f index.sci

genlib:
	cd src/utils && scilab -e "genlib('cgh_utils')"

clean:
	rm -rf temp/ build/ *.log *.bin *.o *.obj

.PHONY: run genlib clean