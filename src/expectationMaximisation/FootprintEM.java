package expectationMaximisation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import main.MIBAT;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import objects.SAMRecordReduced;
import objects.Transcript;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import samTools.SAMReader;
import utils.IO_utils;
import footprintAlignments.FootprintFrameCalculator;

public class FootprintEM {
	public boolean verbose = false;
	private boolean _printProgress = true;
	private void suppressProgressBar(){ _printProgress = false; }
	public boolean printProgressBar(){ return _printProgress; }

	private String _outputPrefix;
	public FootprintEM(String outputPrefix){
		_outputPrefix = outputPrefix;
	}

	//public String recordSummary(SAMRecordReduced record){
	//	return record.getReadName()+"\t"+record.getReferenceName()+"\t"+record.getIntegerAttribute("NH")+"\t"+record.getIntegerAttribute("HI");
	//}


	public void readSAM(File input, boolean restrictToCDS) throws Exception{
		IO_utils.printLineErr("Reading alignments: "+input.getAbsolutePath());
		
		SAMReader engine = new SAMReader(input);

		if(verbose)
			System.err.println("BAM sorted by readID: "+engine.isSortedByReadID());

		if( ! engine.isSortedByReadID()){
			System.err.println("BAM must be sorted by readID ('queryname')\nThis BAM is sorted by: "+engine.getSortOrder().toString());
			System.exit(0);
		}

		ArrayList<SAMRecordReduced> thisRead;
		// Loop through all reads
		while((thisRead = engine.getAlignmentsForNextRead()) != null){

			String thisTranscriptID = "";

			// if this read maps to a single genic locus:
			if(SAMReader.isReadUniquelyMappedToOneGene(thisRead, transcriptID_2_geneID)){

				ArrayList<String> tmp_transcripts = new ArrayList<String>();
				String thisGeneID = "";
				String thisReadID = "";

				// loop again through the alignments and add to the gene list
				for(SAMRecordReduced thisAlignment: thisRead){
					
					// get the transcriptID for this alignment:
					thisTranscriptID = thisAlignment.getReferenceName();

					// if this transcript has a gene in the annotation:
					if(transcriptID_2_geneID.containsKey(thisTranscriptID)){
						thisGeneID = transcriptID_2_geneID.get(thisTranscriptID);
						thisReadID = thisAlignment.getReadName();
						
						Transcript thisTranscript = _genes.get(thisGeneID).getTranscript(thisTranscriptID);
						String transcriptBiotype = thisTranscript.getTranscriptBiotype();

						if(restrictToCDS  &&  _genes.get(thisGeneID).hasCDS()){
							
							double readMid = thisAlignment.getAlignmentStart()+((thisAlignment.getAlignmentEnd()-thisAlignment.getAlignmentStart()+0.0)/2.0);
							
							// If this gene has a CDS, and we only want CDS reads, only allow alignments to coding transcripts
							if(thisTranscript.hasCDS()  &&  transcriptBiotype.equals("protein_coding")  &&  thisTranscript.isCoordInCDS(readMid)){
							//if(thisTranscript.hasCDS()  &&  transcriptBiotype.equals("protein_coding")  &&  thisTranscript.isCoordInCDS(thisAlignment.getAlignmentStart())  &&  thisTranscript.isCoordInCDS(thisAlignment.getAlignmentEnd())){
								tmp_transcripts.add(thisTranscriptID);

								//if(_genes.get(thisGeneID).getgeneSymbol().equals(_printVerboseForGeneName))
								//	System.out.println("has CDS = true,  in CDS = true");

							}//else if(_genes.get(thisGeneID).getgeneSymbol().equals(_printVerboseForGeneName))
							//System.out.println("has CDS = true,  in CDS = FALSE");
						}else{
							// if this gene has no annotated CDS, add all reads 
							tmp_transcripts.add(thisTranscriptID);
							//if(_genes.get(thisGeneID).getgeneSymbol().equals(_printVerboseForGeneName))
							//System.out.println("has CDS = FALSE,  in CDS = FALSE");
						}

					}else{
						System.err.println("ERROR, something went really wrong: transcript \'"+thisTranscriptID+"\' does not appear in the annotation.");
					}
				}

				if(tmp_transcripts.size() > 0){ // read aligns to at least one transcript (may not happen for CDS only filtering)
					// add read to gene:
					this._genes.get(thisGeneID).addRead(thisReadID, tmp_transcripts);
				}

			}else{
				// this is a multi-mapped read.  Ignore.
			}

		}
		
		engine.close();
		IO_utils.printLineErr("Done.");
	}





	/**
	 * 
	 * @param input
	 * @throws Exception
	 */
	public void readGTF(File input) throws Exception{
		IO_utils.printLineErr("Reading GTF: "+input.getAbsolutePath());
		BufferedReader in = new BufferedReader(new FileReader(input));
		String line = "";

		while((line=in.readLine()) != null){
			if(!line.startsWith("##")){
				parseLine(line);
			}
		}
		in.close();
		IO_utils.printLineErr("Done.");
	}

	private HashMap<String, EM_Gene> _genes = new HashMap<String, EM_Gene>();
	private HashMap<String, String> transcriptID_2_geneID = new HashMap<String, String>();


	private void parseLine(String line){
		String[] bits = line.split(" |\t");;

		if(bits[2].equals("exon")  ||  bits[2].equals("CDS")){
			//if(bits[2].equals("CDS")){

			String chromosome = bits[0];
			String strand = bits[6];


			String geneID = trimAttribute(bits[9].trim());
			String transcriptID = trimAttribute(bits[11].trim());
			String geneSymbol = trimAttribute(bits[17].trim());
			String transcriptBiotype = trimAttribute(bits[19].trim());
			String transcriptSymbol = trimAttribute(bits[23].trim());


			// Create gene
			if(!_genes.containsKey(geneID))
				_genes.put(geneID, new EM_Gene(geneID, geneSymbol));

			// Add transcript to gene
			if(!_genes.get(geneID).containsTranscript(transcriptID))
				_genes.get(geneID).addTranscript(transcriptID, chromosome, strand, transcriptBiotype, geneID, transcriptSymbol);

			// Add gene-transcript map
			if(!transcriptID_2_geneID.containsKey(transcriptID))
				transcriptID_2_geneID.put(transcriptID, geneID);

			// Set transcript biotype
			//genes.get(geneID).getTranscript(transcriptID).setBiotype(trimAttribute(bits[19].trim()));



			// Add feature length to the transcript/CDS length
			//int length = Integer.valueOf(bits[4]).intValue()-Integer.valueOf(bits[3]).intValue();
			if(bits[2].equals("CDS")){
				//genes.get(geneID).getTranscript(transcriptID).addCodingExonLength(length);
				_genes.get(geneID).getTranscript(transcriptID).addCDS(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
				_genes.get(geneID).setProteinCoding(true);
			}else{
				//genes.get(geneID).getTranscript(transcriptID).addExonLength(length);
				_genes.get(geneID).getTranscript(transcriptID).addExon(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
			}
		}
	}
	public static String trimAttribute(String in){
		return in.replaceAll("^\"|\";$", "");
	}




	/**
	 * 
	 * @param inputFile
	 * @throws Exception
	 */
	public void readTranscriptExpressions(File inputFile) throws IOException{
		IO_utils.printLineErr("Reading eXpress expression data: "+inputFile.getAbsolutePath());
		BufferedReader in = new BufferedReader(new FileReader(inputFile));
		String line = "";
		String[] bits;

		// burn the first line
		in.readLine();

		while((line=in.readLine()) != null){
			bits = line.split("\t");

			//transcriptQuants.put(bits[1].trim(), Double.valueOf(bits[14].trim()).doubleValue());
			String transcriptID = bits[1].trim();
			if(this.transcriptID_2_geneID.containsKey(transcriptID))
				this._genes.get(this.transcriptID_2_geneID.get(transcriptID)).getTranscript(transcriptID).setTPM(Double.valueOf(bits[14].trim()).doubleValue());;
		}
		in.close();

		IO_utils.printLineErr("Done.");
	}


	/**
	 * 
	 */
	public void initialisePriors(){
		//System.out.println(this.getTime()+" Initialising priors...");

		Iterator<String> iterator = _genes.keySet().iterator();
		String geneID = "";

		while(iterator.hasNext()){
			geneID = iterator.next();
			_genes.get(geneID).setPriors();
		}

		//System.out.println(this.getTime()+" Done.");
	}

	private String _printVerboseForGeneName = "";

	private HashMap<String, String> _readMap = new HashMap<String, String>();
	public void doEM(int maxIterations, double converganceDistance, boolean outputAll, boolean restrictToCDS) throws IOException{
		//System.err.println(getTime()+" Running EM (maxIterations: "+maxIterations+"), writing to: "+this._outputPrefix+".exprs");
		IO_utils.printLineErr("Running EM (maxIterations: "+maxIterations+"), writing to: "+this._outputPrefix+".exprs");
		Iterator<String> iterator = _genes.keySet().iterator();
		int totalGenes = _genes.size();
		int count = 0;
		double percent = 0.0;
		String geneID = "";
		EM_Core_optimised emCore;

		PrintWriter out = new PrintWriter(new FileWriter(this._outputPrefix+".exprs"));
		//System.out.println("geneID\ttranscriptID\thasFlatPrior\tprior\tflag\tnIterations\teffectiveReadCount\treadDensity\ttranscriptFractionAfterFootprints");
		out.println("geneID\ttranscriptID\ttranscriptLength_exons\ttranscriptLength_CDS\tnFootprintsMappedToThisGene\tnReadsMappedToThisTranscript\thasFlatPrior\tprior\tflag\tnIterations\teffectiveReadCount\treadDensity_transcript\treadDensity_CDS\ttranscriptFractionAfterFootprints\tgeneSymbol\ttranscriptName\ttranscriptBiotype");

		while(iterator.hasNext()){
			geneID = iterator.next();
			emCore = new EM_Core_optimised(_genes.get(geneID), restrictToCDS);

			if(_genes.get(geneID).getgeneSymbol().equals(_printVerboseForGeneName))
				emCore.setVerbose(true);

			//if(_gene2inframeRead2isoform.containsKey(geneID))
			//	emCore.filterRead2TranscriptAlignments(_gene2inframeRead2isoform.get(geneID));

			EM_Result res = emCore.runEM(maxIterations, converganceDistance);//, out);
			_readMap = emCore.sampleReadAssignments(_readMap);

			// Format and print EM results
			if(outputAll  ||  res.getGene().getReadCounts().size() > 0){
				for(int i=0;i<res.getGene().getTranscriptIDs().size();i++){
					String thisTranscriptID = res.getGene().getTranscriptIDs().get(i);
					double tmp_effectiveCount = res.getFinalEffectiveReadCount().get(thisTranscriptID).doubleValue()+0.0;
					//double tmp_coverage = tmp_effectiveCount / res.getGene().getLengths().get(thisTranscriptID);
					double tmp_finalLikelihood = res.getFinalTranscriptLikelihoods().get(thisTranscriptID).doubleValue();
					int exonlength = res.getGene().getTranscript(thisTranscriptID).getTotalExonLength();
					int cdslength = res.getGene().getTranscript(thisTranscriptID).getTotalCodingExonLength();
					double density_exon = tmp_effectiveCount/(exonlength+0.0);
					double density_cds = tmp_effectiveCount/(cdslength+0.0);
					if(tmp_effectiveCount == 0.0)
						density_exon = density_cds = 0.0;
					else if(cdslength == 0)
						density_cds = 0.0;
					out.printf(geneID+"\t"+
							thisTranscriptID+"\t"+
							exonlength+"\t"+
							cdslength+"\t"+
							res.getGene().getReadCounts().size()+"\t"+
							res.getGene().getTranscriptReadCounts(thisTranscriptID)+"\t"+
							res.getGene().hasFlatPrior()+"\t%e\t"+
							res.getFlag()+"\t"+
							res.getNIterations()+"\t%f\t%e\t%e\t%e\t"+
							res.getGene().getgeneSymbol()+"\t"+
							res.getGene().getTranscript(thisTranscriptID).getTranscriptSymbol()+"\t"+
							res.getGene().getTranscript(thisTranscriptID).getTranscriptBiotype()+"\n",
							res.getGene().getPriors().get(thisTranscriptID), tmp_effectiveCount, density_exon, density_cds, tmp_finalLikelihood);
				}
			}			
			count++;

			if(this.printProgressBar()){
				if(Math.round((count*100.0)/totalGenes) > percent)
					percent = (int)Math.round((count*100.0)/totalGenes);
				IO_utils.printProgressBar(percent);
			}

		}

		if(this.printProgressBar())
			IO_utils.printErr("\n");

		out.flush();
		out.close();
		//System.err.println(this.getTime()+" Done.");
		IO_utils.printLineErr("Done.");
	}


	/**
	 * 
	 * @param readMap
	 * @throws IOException
	 */
	public void outputSampleReadAssignments(File inputSAM) throws IOException{
		//PrintWriter out = new PrintWriter(new FileWriter(this._outputPrefix+".bam"));
		IO_utils.printLineErr("Writing sampled read alignments to: "+this._outputPrefix+".sample.bam");

		SAMFileReader inputSam = new SAMFileReader(inputSAM);
		SAMFileWriterFactory fact=new SAMFileWriterFactory();
		fact.setCreateIndex(true);
		SAMFileWriter outputSam = fact.makeBAMWriter(inputSam.getFileHeader(), false, new File(this._outputPrefix+".sample.bam"));

		SAMRecord thisRecord;
		SAMRecordIterator it = inputSam.iterator();
		String thisTranscriptID,thisReadID; 
		while(it.hasNext()){
			thisRecord = it.next();
			thisTranscriptID = thisRecord.getReferenceName();
			thisReadID = thisRecord.getReadName();

			if(_readMap.containsKey(thisReadID)){
				if(thisTranscriptID.equals(_readMap.get(thisReadID))){
					outputSam.addAlignment(thisRecord);
					_readMap.remove(thisReadID);
				}
			}
		}
		outputSam.close();
		inputSam.close();
		IO_utils.printLineErr("Done.");
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("annotationPath").hasArg().withDescription("Path to the GTF file containing the transcript/gene relationship").create(MIBAT.OPT_PATH_ANNOTATION));
		options.addOption(OptionBuilder.withArgName("results.xprs").hasArg().withDescription("[optional] Text file containing RNA-seq transcript expression quantifications from eXpress").create("e"));
		options.addOption(OptionBuilder.withArgName("footprint alignments").hasArg().withDescription("SAM/BAM file containing alignments against the transcriptome").create("f"));
		options.addOption(OptionBuilder.withLongOpt("outputPrefix").withArgName("outputPath").hasArg().withDescription("Output prefix for results").create("o"));
		options.addOption(OptionBuilder.withLongOpt("sample").withDescription("Output sampling of read alignments based on EM transcript likelihoods").create("s"));
		options.addOption(OptionBuilder.withArgName("maxIterations_EM").hasArg().withDescription("Maximum number of EM iterations [default: 500]").create("N"));
		options.addOption(OptionBuilder.withArgName("double").hasArg().withDescription("Criteria for EM convergence [default: 0.0001]").create("c"));
		options.addOption(OptionBuilder.withDescription("Write all genes/transcripts, even those with no observed footprints").create("v"));
		options.addOption(OptionBuilder.withDescription("Count only footprints in the CDS (where available)").create("cds"));
		options.addOption(OptionBuilder.withDescription("Include frame information for the footprints in the EM [automatically adds --cds]").create("useFrame"));
		options.addOption(OptionBuilder.withDescription("Keep intermediate files").create("DEV"));
		options.addOption(OptionBuilder.withDescription("Do not print progress to stderr").create("noprog"));
		return options;
	}


	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {

		/*args = new String[]{"IsoformEM_Footprints", 
						//"-f", "/Users/robk/Desktop/EM_TEST/TEST_NEW/D1-Footprint_1x50_gencode.v18.annotation_mappedReads.sorted.bam",
						"-f", "/Users/robk/Desktop/EM_TEST/TEST_NEW/ALL-Footprint_1x50_gencode.v18.annotation_mappedReads.sorted.bam",
						"-a", "/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v18.annotation.filtered.gtf",
						"-e", "/Users/robk/Box Sync/Work/HEK293_RNAseq/STAR_eXpress/A1-Total_2x75/results.xprs",
						"-o", "/Users/robk/Desktop/EM_TEST/TEST_NEW",
						"-DEV"};
		 */
		/*
		boolean useFrame = true;

		String frameFlag = "";
		if(useFrame)
			frameFlag = "--useFrame";

		args = new String[]{"IsoformEM_Footprints",
				"-a","/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v21.annotation_noSelenocysteine.gtf",
				"-f","/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Footprinting/FP_original_1_Genome_Aligned.toTranscriptome.sorted.bam", 
				"-N","1000",
				"--sample",
				"--cds", 
				frameFlag,
				"--noprog",
				"-o","/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Footprinting/TESTTEST_frame-"+useFrame};
		 */

		CommandLine cmdArgs = MIBAT.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption(MIBAT.OPT_PATH_ANNOTATION) && cmdArgs.hasOption("f") && cmdArgs.hasOption("o")){
			System.err.println();

			// initialise
			FootprintEM engine = new FootprintEM(cmdArgs.getOptionValue("o"));

			// for testing!
			//engine._printVerboseForGeneName = "COG8";
			//engine._printVerboseForGeneName = "POLDIP3";

			// suppress printing of the progress bar
			if(cmdArgs.hasOption("noprog"))
				engine.suppressProgressBar();

			// should we restrict footprints to the CDS of coding genes
			boolean restrictToCDS = cmdArgs.hasOption("cds");

			// status update
			IO_utils.printLineErr("Include frame information in EM: "+cmdArgs.hasOption("useFrame")); 
			// should we include footprint frame information in the EM
			if(cmdArgs.hasOption("useFrame")){
				engine.performFrameCalc(cmdArgs.getOptionValue("f"), cmdArgs.getOptionValue(MIBAT.OPT_PATH_ANNOTATION));
				restrictToCDS = true;
			}

			// status update			
			IO_utils.printLineErr("Restrict to CDS aligned reads: "+restrictToCDS);

			// Read the transcripts in the GTF
			engine.readGTF(new File(cmdArgs.getOptionValue(MIBAT.OPT_PATH_ANNOTATION)));

			// Read the eXpress transcript quantifications or choose flat priors
			if(cmdArgs.hasOption("e")){
				engine.readTranscriptExpressions(new File(cmdArgs.getOptionValue("e")));
			}

			// Convert the transcript expressions to priors and output
			engine.initialisePriors();

			// Parse the footprint alignments
			engine.readSAM(new File(cmdArgs.getOptionValue("f")), restrictToCDS);

			// modify the alignment weights for the EM using frame information
			if(cmdArgs.hasOption("useFrame"))
				engine.readFootprintFrames();
			else
				engine.resetAlignmentWeights();

			// set the number of iterations to cap the EM
			int maxIterations = 500;
			if(cmdArgs.hasOption("N"))
				maxIterations = Integer.valueOf(cmdArgs.getOptionValue("N")).intValue();

			// set the convergence criteria for the EM
			double convergenceDistance = 1.0/10000.0;
			if(cmdArgs.hasOption("c"))
				convergenceDistance = Double.valueOf(cmdArgs.getOptionValue("c")).doubleValue();

			boolean outputAll = false;
			if(cmdArgs.hasOption("v"))
				outputAll = true;

			// Run the EM!
			//engine.runEM("perl -I "+tmp_path+" "+tmp_path+"/main.foot.print.pl "+maxIterations+" "+output_priors+" "+output_alignments+" "+output_path);
			engine.doEM(maxIterations, convergenceDistance, outputAll, restrictToCDS);

			// Sample the reads?
			if(cmdArgs.hasOption("s")){
				engine.outputSampleReadAssignments(new File(cmdArgs.getOptionValue("f")));
			}

			// Clean up
			if(cmdArgs.hasOption("DEV")){
				//engine.removeFile(new File(tmp_path));
				engine.verbose = true;
			}

			IO_utils.printLineErr("All Done!");

		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(MIBAT.THUNDER_EXE_COMMAND+" IsoformEM_Footprints", getCmdLineOptions());
			System.err.println();
		}


	}

	/** stores read/transcript alignment only for those reads that are in frame with the CDS*/
	//HashMap<String, HashMap<String, HashMap<String, Integer>>> _gene2isoform2read2inframe = new HashMap<String, HashMap<String, HashMap<String, Integer>>>();
	//HashMap<String, HashMap<String, HashSet<String>>> _gene2isoform2inframeRead = new HashMap<String, HashMap<String, HashSet<String>>>();
	HashMap<String, HashMap<String, HashSet<String>>> _gene2inframeRead2isoform = new HashMap<String, HashMap<String, HashSet<String>>>();


	HashMap<Integer, HashMap<Double, Double>> _readLength2offset2weight = new HashMap<Integer, HashMap<Double, Double>>();


	private void resetAlignmentWeights(){
		Iterator<String> it = _genes.keySet().iterator();
		while(it.hasNext())
			_genes.get(it.next()).resetAlignmentWeights(1.0);
	}

	/**
	 * Function to call the frame analysis module 
	 * @param readsPath
	 * @param annotationPath
	 * @param outputPrefix
	 * @throws Exception
	 */
	private void performFrameCalc(String readsPath, String annotationPath) throws Exception{
		FootprintFrameCalculator engine = new FootprintFrameCalculator(annotationPath, 1);
		if(!_printProgress)
			engine.suppressProgressBar();
		engine.doAll(new File(readsPath), 1000000000, _outputPrefix, false, true);

		// Get global read-length to offset
		_readLength2offset2weight = engine.getReadLengths2offsets2weights();
		IO_utils.printLineErr("# read lengths:"+_readLength2offset2weight.size());
	}


	/**
	 * Reads the results of the frame calculation back to the HashMap above
	 * @param outputPrefix
	 * @throws Exception
	 */
	private void readFootprintFrames() throws Exception{
		int countOfInvalidReadFrames = 0;
		//
		// Read individual read frame offsets
		//
		BufferedReader in = new BufferedReader(new FileReader(new File(_outputPrefix+"_Read2TranscriptFrame.txt")));
		//String header = in.readLine(); //burn header line
		in.readLine(); //burn header line
		String line = "";
		String[] bits;
		while((line=in.readLine())!=null){
			bits = line.split("\t");

			String geneID = bits[0].trim();
			String transcriptID = bits[1].trim();
			String readID = bits[2].trim();
			int readLength = Integer.valueOf(bits[3].trim()).intValue();
			boolean readSense = Boolean.valueOf(bits[4].trim()).booleanValue();
			double frameOffset = Double.valueOf(bits[5].trim()).doubleValue();
			String inCDS = bits[6].trim();

			if(_genes.get(geneID).getgeneSymbol().equals(_printVerboseForGeneName))
				//System.out.print("\n"+header+"\n"+line);
				System.out.print(line);
				
			if(inCDS.equals("true")  &&  readSense == true  &&  frameOffset >= 0.0  &&  _genes.containsKey(geneID)){
				if(_genes.get(geneID).containsTranscript(transcriptID)){
					try{
						_genes.get(geneID).addWeightToReadAlignment(readID, transcriptID, _readLength2offset2weight.get(readLength).get(frameOffset));
						
						if(_genes.get(geneID).getgeneSymbol().equals(_printVerboseForGeneName))
							System.out.print("\t--> using");
						
					}catch(NullPointerException e){
						//System.err.println(geneID+"\t"+transcriptID+"\t"+readID+"\t"+readLength+"\t"+frameOffset+"\t"+_readLength2offset2weight.get(readLength).containsKey(frameOffset));
						countOfInvalidReadFrames ++;
					}
				}
			}
			
			if(_genes.get(geneID).getgeneSymbol().equals(_printVerboseForGeneName))
				System.out.println();
		}
		in.close();
		IO_utils.printLineErr("Invalid frame offset detected in "+countOfInvalidReadFrames+" alignmnments");
	}






	public void makeOutputDirectories(String output_path){
		File tmp = new File(output_path);
		if(!tmp.exists())
			tmp.mkdir();
		//tmp = new File(output_path+"/tmp");
		//if(!tmp.exists())
		//	tmp.mkdir();
	}

	public void removeFile(File file){
		if(file.isDirectory()){
			File[] children = file.listFiles();
			for(int i=0;i<children.length;i++)
				removeFile(children[i]);
		}
		file.delete();
		//System.out.println("Deleted: "+file.getAbsolutePath());
	}

}



