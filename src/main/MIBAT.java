package main;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import proteome.ProcessPepDigest;
import expectationMaximisation.FootprintEM;
import expectationMaximisation.ProteomicEM;
import footprintAlignments.FootprintFrameCalculator;

public class MIBAT {

	public static final String VERSION = "0.6.1";
	
	public static final String OPT_PATH_DB_ANNOTATION = "A";
	public static final String OPT_PATH_DB_SAMPLE = "S";
	public static final String OPT_DB_TABLE_NAME = "t";
	public static final String OPT_PATH_INPUT = "i";
	public static final String OPT_PATH_OUTPUT = "o";
	public static final String OPT_PATH_SPECTRA_TANDEM = "i";
	public static final String OPT_PATH_SPECTRA_MZXML = "i";
	public static final String OPT_PATH_ANNOTATION = "a";
	public static final String OPT_FORMAT_GENCODE = "g";
	public static final String OPT_FORMAT_CUFFMERGE = "c";
	public static final String OPT_FORMAT_SWISSPROT = "s";
	public static final String OPT_FORMAT_TRINITY = "n";
	public static final String OPT_CHOICE_DB_FORCE_REFRESH = "F";

	/**
	 * Parse the command line arguments
	 * @param args
	 * @return
	 * @throws ParseException
	 */
	public static CommandLine parseArgs(String[] args, Options options) throws ParseException{
		CommandLineParser parser = new BasicParser();
		return parser.parse(options, args);	
	}



	public static void main(String[] args) throws Exception {
		String main = "NULL";
		if(args.length > 0){
			main = args[0].trim().toLowerCase(); 
		}

		if(main.equals("isoformem_footprints")){
			FootprintEM.main(args);
		}else if(main.equals("isoformem_proteomics")){
			ProteomicEM.main(args);
		}else if(main.equals("footprintframeanaysis")){
			FootprintFrameCalculator.main(args);
		}else if(main.equals("peptidedigest")){
			ProcessPepDigest.main(args);
		}
		
		
		
		else{
			System.out.println("EMpire version "+VERSION);
			System.out.println("");
			
			System.out.println("Usage:\t "+MIBAT.THUNDER_EXE_COMMAND+" <Command>");
			System.out.println("");

			System.out.println("Command: IsoformEM_Footprints  | Infer most likely transcripts from ribosome footprint alignments");
			System.out.println("         IsoformEM_Proteomics  | Infer most likely isoforms from MS/MS spectra mapping");
			System.out.println("         FootprintFrameAnaysis | Analyse ribosome footprint alignments in terms of fidelity to annotated coding frames");
			System.out.println("         PeptideDigest         | Enzymatically digest a protein reference");
			System.out.println();
		}


	}

	public static final String THUNDER_EXE_COMMAND = "java -Xmx2G -jar EMpire.jar";

	
	
}
