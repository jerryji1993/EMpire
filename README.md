# EMpire

Simply download and run EMpire.jar for a list of available sub-commands.

Usage:	 java -Xmx2G -jar EMpire.jar \<Command\>

Available \<Command\> options: 

	IsoformEM_Footprints  | Infer most likely transcripts from ribosome footprint alignments
 	IsoformEM_Proteomics  | Infer most likely isoforms from MS/MS spectra mapping
	FootprintFrameAnaysis | Analyse ribosome footprint alignments in terms of fidelity to annotated coding frames
	PeptideDigest         | Enzymatically digest a protein reference


## Depedencies

Relies heavily on general purpose java code in https://github.com/rkitchen/Thunder

Also requires several java libraries:
- commons-io-2.4
- commons-cli-1.2
- commons-lang3-3.3.2
- commons-math3-3.5
- sqlite-jdbc-3.7.2
- sam-1.96
