NCBI_SERVER = ftp.ncbi.nlm.nih.gov
FTP_SERVER = ftp://$(NCBI_SERVER)
RSYNC_SERVER = "rsync://$NCBI_SERVER"

clean:
	rm -f .complete
	find . -mindepth 1 -type d |  xargs rm -fr
