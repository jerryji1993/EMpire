package expectationMaximisation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import main.MIBAT;
import objects.Alignment;
import objects.MS1_scan;
import objects.Peptide;
import objects.Spectra;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import proteome.ProcessPepDigest;
import proteome.ReadXTandem;
import proteome.SAXHandler_mzXML;
import proteome.SpectraAlignmentEngine;
import utils.IO_utils;

public class ProteomicEM {
	public boolean verbose = false;

	private String _outputPrefix;
	public ProteomicEM(String outputPrefix){
		_outputPrefix = outputPrefix;
	}

	private String _printVerboseForGeneName = "";


	/*private HashMap<String, HashMap<Integer,Integer>> _isoformPeptideCounts;
	public void readDigestInfo(File input) throws IOException{
		System.err.println(this.getTime()+" Reading peptide digest info: "+input.getAbsolutePath());
		_isoformPeptideCounts = ProcessPepDigest.processEMBOSSoutput(input);
		System.err.println(this.getTime()+" Done.");
	}*/
	private HashMap<String, Integer> _isoformPeptideCounts = new HashMap<String, Integer>();
	public void readDigestInfo(File input) throws IOException{
		IO_utils.printLineErr("Performing in-silico peptide digest: "+input.getAbsolutePath());
		ArrayList<String> enzymes = new ArrayList<String>();
		enzymes.add(ProcessPepDigest.ENZYME_TRYPSIN);
		enzymes.add(ProcessPepDigest.ENZYME_LYSC);
		//_isoformPeptideCounts = ProcessPepDigest.digestProteinsFromFasta(input, enzymes);

		HashMap<String, Integer> tmpPeptides = ProcessPepDigest.digestProteinsFromFasta(input, enzymes);

		Iterator<String> it = tmpPeptides.keySet().iterator();
		String tmpID;
		String[] isoformFrameGene;
		while(it.hasNext()){
			tmpID = it.next();
			isoformFrameGene = SpectraAlignmentEngine.parseIsoformID(tmpID);
			_isoformPeptideCounts.put(isoformFrameGene[0]+"_"+isoformFrameGene[1], tmpPeptides.get(tmpID));
		}

		IO_utils.printLineErr("Done.");
	}


	private HashMap<String, MS1_scan> _ms1Intensities = new HashMap<String, MS1_scan>();
	/**
	 * 
	 * @param input
	 */
	public void readMZXML(File input){
		IO_utils.printLineErr("Reading spectra intensities (MS1): "+input.getAbsolutePath());

		try {
			SAXParserFactory factory = SAXParserFactory.newInstance();
			SAXParser saxParser = factory.newSAXParser();
			SAXHandler_mzXML handler = new SAXHandler_mzXML();

			saxParser.parse(input, handler);
			_ms1Intensities = handler.getMS1Intensities();
			//System.err.println();

			IO_utils.printLineErr("N spectra: "+_ms1Intensities.size());

		}catch (Exception e) {
			e.printStackTrace();
		}finally{
		}

		/*Iterator<String> it = ms1Intensities.keySet().iterator();
		for(int i=0;i<10;i++){
			System.out.println(it.next());
		}*/

		IO_utils.printLineErr("Done.");
	}



	/**
	 * 
	 * @param input
	 * @throws Exception
	 */
	public void readSpectra_tandem(String tandemResultsPath, double fdrMax, boolean restrictToCDS) throws Exception{
		String[] filePaths = tandemResultsPath.split(",");
		if(filePaths.length >= 1){
			File[] inputFiles = new File[filePaths.length];
			for(int i=0;i<filePaths.length;i++){
				inputFiles[i] = new File(filePaths[i].trim());
			}

			// Read spectra
			//IO_utils.printLineErr("Reading spectra alignments: "+tandemResults.getAbsolutePath());
			SpectraAlignmentEngine spectraAlignments = ReadXTandem.processTandemXML(inputFiles, _outputPrefix, fdrMax, true);

			String thisGeneID = "";
			String thisTranscriptID = "";
			String thisSpectraID = "";
			String thisFrame = "";
			ArrayList<String> tmp_transcripts = new ArrayList<String>();

			Alignment currentAlignment = null;
			Iterator<String> allSpectra = spectraAlignments.getAlignments().keySet().iterator();
			Iterator<Alignment> currentSpectraAlignments;

			ArrayList<String> tmp_genes;

			String[] isoformFrameGene;
			ArrayList<String> unknownIsoformIDs = new ArrayList<String>();
			while(allSpectra.hasNext()){
				currentSpectraAlignments = spectraAlignments.getAlignments().get(allSpectra.next()).iterator();

				int spectraIntensity = 1;
				boolean hasMS1 = false;

				tmp_genes = new ArrayList<String>();
				tmp_transcripts = new ArrayList<String>();

				while(currentSpectraAlignments.hasNext()){
					currentAlignment = currentSpectraAlignments.next();

					isoformFrameGene = SpectraAlignmentEngine.parseIsoformID(currentAlignment.getProtein().getIsoformID());
					thisTranscriptID = isoformFrameGene[0];
					thisGeneID = isoformFrameGene[2];
					thisFrame = isoformFrameGene[1];
					//thisSpectraID = currentAlignment.getSpectra().getAttribute(Spectra.ATTRIBUTE_ID);
					thisSpectraID = currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_SEQ);

					//if(_genes.get(thisGeneID).getgeneSymbol().equals(_printVerboseForGeneName))
					//	System.err.println("\n>readSpectra");

					if(transcriptID_2_geneID.containsKey(thisTranscriptID+"_"+thisFrame)){					
						thisGeneID = transcriptID_2_geneID.get(thisTranscriptID+"_"+thisFrame);

						if(!tmp_genes.contains(thisGeneID))
							tmp_genes.add(thisGeneID);

						//if(_genes.get(thisGeneID).getgeneSymbol().equals(_printVerboseForGeneName))
						//	System.out.print("geneID="+thisGeneID+"\ttranscriptID="+thisTranscriptID+"\tframe="+thisFrame+"\trestrictToCDS="+restrictToCDS+"\tgeneHasCDS="+_genes.get(thisGeneID).hasCDS());


						/*if(thisGeneID.equals("ENSG00000181652.16")){
						System.out.println("Gene has CDS: "+this.genes.get(thisGeneID).hasCDS());
						System.out.println(thisTranscriptID+" biotype: "+this.genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).getBiotype());
						System.out.println(thisTranscriptID+" has CDS: "+this.genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).hasCDS());
					}*/

						if(restrictToCDS  &&  _genes.get(thisGeneID).hasCDS()){
							// If this gene has a CDS, and we only want CDS peptides, only allow alignments to coding transcripts
							int peptideStart = Integer.valueOf(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_START)).intValue();
							int peptideEnd = Integer.valueOf(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_END)).intValue();
							if(_genes.get(thisGeneID).getgeneSymbol().equals(_printVerboseForGeneName))
								System.err.print("Reading spectra: "+
										""+thisSpectraID+"\tgeneID="+thisGeneID+"\ttranscriptID="+thisTranscriptID+"_"+thisFrame+
										"\ttranscriptHasCDS="+_genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).hasCDS()+
										"\tbiotype="+_genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).getTranscriptBiotype()+
										"\tpeptideStart="+peptideStart+"\tpeptideEnd="+peptideEnd+
										"\ttranscript_Exonstart="+_genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).getTranscriptStart()+
										"\ttranscript_CDSstart="+_genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).getCDSStart()+
										"\ttranscript_CDSstop="+_genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).getCDSStop());

							if(_genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).hasCDS()  &&
									_genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).getTranscriptBiotype().equals("protein_coding")){//  &&
								//peptideStart >= _genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).getCDSStart()  &&
								//peptideEnd <= _genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).getCDSStop()){

								if(_genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).isCoordInCDS(peptideStart)  &&  
										_genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).isCoordInCDS(peptideEnd)){

									tmp_transcripts.add(thisTranscriptID+"_"+thisFrame);

									if(_genes.get(thisGeneID).getgeneSymbol().equals(_printVerboseForGeneName))
										System.err.print("\t--> ADDED");
								}
							}
						}else{
							// if this gene has no annotated CDS, add all peptides
							tmp_transcripts.add(thisTranscriptID+"_"+thisFrame);
						}

						if(_genes.get(thisGeneID).getgeneSymbol().equals(_printVerboseForGeneName))
							System.err.println();


					}else{
						// this transcript is not in the GTF! 
						if(!unknownIsoformIDs.contains(currentAlignment.getProtein().getIsoformID())){	
							unknownIsoformIDs.add(currentAlignment.getProtein().getIsoformID());
							System.err.println("WARNING: Unable to find annotation information for "+currentAlignment.getProtein().getIsoformID());
						}
					}

					if(!hasMS1){
						if(_ms1Intensities.containsKey(currentAlignment.getSpectra().getAttribute(Spectra.ATTRIBUTE_ID))){
							hasMS1 = true;
							spectraIntensity = (int)Math.round(Double.valueOf(_ms1Intensities.get(currentAlignment.getSpectra().getAttribute(Spectra.ATTRIBUTE_ID)).getPrecursorAttribute(MS1_scan.ATTRIBUTE_PRECURSOR_INTENSITY)).doubleValue());
							//System.out.println(currentAlignment.getSpectra().getAttribute(Spectra.ATTRIBUTE_ID)+": "+spectraIntensity);
						}
					}
				}

				// add peptides to gene if there is only one gene:
				if(tmp_genes.size() == 1  &&  tmp_transcripts.size() > 0){
					if(_genes.get(thisGeneID).getgeneSymbol().equals(_printVerboseForGeneName)){
						System.err.print("Adding "+thisSpectraID+"\tto "+thisGeneID+"\ttranscripts:\t");
						for(String transcript: tmp_transcripts)
							System.err.print(transcript+"\t");

						System.err.println();
						//_genes.get(thisGeneID).printInfo();
					}

					if(hasMS1)
						this._genes.get(thisGeneID).addRead(thisSpectraID, tmp_transcripts, spectraIntensity);
					else
						this._genes.get(thisGeneID).addRead(thisSpectraID, tmp_transcripts);

					this._genes.get(thisGeneID).resetAlignmentWeights(1.0);
				}

			}
			IO_utils.printLineErr("Done.");
		}
		else{
			IO_utils.printLineErr("ERROR: no input files specified");
		}
	}



	/**
	 * 
	 * @param maxquantPeptidesPath
	 * @param sampleID
	 * @return
	 * @throws IOException
	 */
	private boolean readPeptides_maxquant(String maxquantPeptidesPath, String sampleID) throws IOException{
		BufferedReader in = new BufferedReader(new FileReader(maxquantPeptidesPath));

		// look at the header - do we have the sample prefix we're expecting?
		int colIndex = -1;
		String detectedSampleIDs = "";
		String line = in.readLine();
		if(line.startsWith("Sequence")){
			String[] linebits = line.split("\t");
			for(int i=0;i<linebits.length; i++){
				if(linebits[i].startsWith("Intensity")){
					// if this is one of the columns we might use, check for the correct sample ID
					detectedSampleIDs = detectedSampleIDs.concat("\'"+linebits[i].substring(9).replace(".", "")+"\' ");		
					if(linebits[i].endsWith(sampleID)){
						//System.out.println(linebits[i]);
						colIndex = i;
					}
				}
			}
		}

		// If we have detected the correct sample, get the mapping results
		if(colIndex >= 0){
			while((line=in.readLine())!=null){
				//System.out.println(line);
				String[] linebits = line.split("\t");
				String peptideSeq = linebits[0];
				String[] transcriptIDs_raw = linebits[34].split(";");

				// check if any of the transcripts are contaminants - if so, skip this peptide
				boolean usePeptide = false;
				for(int i=0;i<transcriptIDs_raw.length;i++){
					if(transcriptIDs_raw[i].startsWith("ENS")){
						usePeptide = true;
					}else if(transcriptIDs_raw[i].startsWith("CON")){
						usePeptide = false;
						break;
					}
				}

				if(usePeptide){
					try{
						long intensity = Double.valueOf(linebits[colIndex]).longValue();
						//long intensity = Long.parseLong(linebits[colIndex]);
						//int intensity = Integer.valueOf(linebits[colIndex]).intValue();

						if(intensity > 0){
							String[] transcriptIDs = new String[transcriptIDs_raw.length];
							String[] transcriptFrames = new String[transcriptIDs_raw.length];

							// split the transcripts to which this peptide maps
							for(int i=0;i<transcriptIDs_raw.length;i++){
								String[] tmp = transcriptIDs_raw[i].split("_");
								transcriptIDs[i] = tmp[0];
								transcriptFrames[i] = tmp[1];
								if(transcriptFrames[i].length() == 0)
									System.err.println(transcriptIDs_raw[i]);
							}
							//System.out.print(peptideSeq+"\t"+intensity+"\t");

							// loop through all transcript alignments 
							ArrayList<String> tmp_genes = new ArrayList<String>();
							ArrayList<String> tmp_transcripts = new ArrayList<String>();
							for(int i=0;i<transcriptIDs.length;i++){
								String thisTranscriptID = transcriptIDs[i];
								
								int thisFrame = Integer.valueOf(transcriptFrames[i]).intValue();
								//System.out.print(transcriptIDs[i]+":"+transcriptFrames[i]+"\t");
								if(transcriptID_2_geneID.containsKey(thisTranscriptID+"_"+thisFrame)){					
									String thisGeneID = transcriptID_2_geneID.get(thisTranscriptID+"_"+thisFrame);

									if(!tmp_genes.contains(thisGeneID))
										tmp_genes.add(thisGeneID);

									tmp_transcripts.add(thisTranscriptID+"_"+thisFrame);

									if(_genes.get(thisGeneID).getgeneSymbol().equals(_printVerboseForGeneName))
										System.err.print("\t--> ADDED");
								}

							}
							//System.out.println();

							// add peptides to gene if there is only one gene:
							//System.out.println(tmp_genes.size());
							if(tmp_genes.size() == 1  &&  tmp_transcripts.size() > 0){
								String thisGeneID = tmp_genes.get(0);
								//						if(_genes.get(thisGeneID).getgeneSymbol().equals(_printVerboseForGeneName)){
								//							System.err.print("Adding "+thisSpectraID+"\tto "+thisGeneID+"\ttranscripts:\t");
								//							for(String transcript: tmp_transcripts)
								//								System.err.print(transcript+"\t");
								//
								//							System.err.println();
								//							//_genes.get(thisGeneID).printInfo();
								//						}



								//this._genes.get(thisGeneID).addRead(peptideSeq, tmp_transcripts, 1);
								this._genes.get(thisGeneID).addRead(peptideSeq, tmp_transcripts, intensity);

								this._genes.get(thisGeneID).resetAlignmentWeights(1.0);
							}

						}
					}catch(Exception e){
						IO_utils.printLineErr("ERROR parsing line: "+line);
						e.printStackTrace(System.err);
					}
				}
			}
		}
		in.close();

		// if we don't have a column index, 
		if(colIndex < 0){
			if(sampleID == null){
				IO_utils.printLineErr("No sample ID specified, detected the following in the file - please choose one of the following and re-run EMpire with the --sampleID {sampleID} flag.");
			}else{
				IO_utils.printLineErr("Unable to find specified sampleID, detected the following in the file - please choose one of the following and re-run EMpire.");
			}
			IO_utils.printLineErr(detectedSampleIDs);
			return false;
		}else
			return true;
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

		for(String thisGene: _genes.keySet()){
			for(String thisTranscript: _genes.get(thisGene).getTranscriptIDs()){
				if(_isoformPeptideCounts.containsKey(thisTranscript)){
					_genes.get(thisGene).getTranscript(thisTranscript).setExonLength(_isoformPeptideCounts.get(thisTranscript));
					_genes.get(thisGene).getTranscript(thisTranscript).setCodingExonLength(_isoformPeptideCounts.get(thisTranscript));
				}else{
					_genes.get(thisGene).getTranscript(thisTranscript).setExonLength(0);
					_genes.get(thisGene).getTranscript(thisTranscript).setCodingExonLength(0);
				}
			}
		}

		IO_utils.printLineErr("Done.");
	}

	private HashMap<String, EM_Gene> _genes = new HashMap<String, EM_Gene>();
	private HashMap<String, String> transcriptID_2_geneID = new HashMap<String, String>();




	/*
	 * Read GTF entries to memory - IMPORTANT: read all transcripts in THREE FRAMES!
	 */
	private void parseLine(String line){
		String[] bits = line.split(" |\t");;

		String featureType = bits[2].trim(); 
		if(featureType.equals("exon")  ||  featureType.equals("CDS")){

			String chromosome = bits[0];
			String strand = bits[6];

			String geneID = trimAttribute(bits[9].trim());
			String transcriptID = trimAttribute(bits[11].trim());
			String geneSymbol = trimAttribute(bits[17].trim());
			String transcriptBiotype = trimAttribute(bits[19].trim());
			String transcriptSymbol = trimAttribute(bits[23].trim());


			//			if(geneID.equals("ENSG00000181652.16")){
			//				System.out.println(featureType);
			//			}

			// Create gene
			if(!_genes.containsKey(geneID)){
				_genes.put(geneID, new EM_Gene(geneID, geneSymbol));
			}


			/*// Placeholder variables
			int tmp_peptideCount_1 = 0;
			int tmp_peptideCount_2 = 0;
			int tmp_peptideCount_3 = 0;
			 */

			// Add transcript in all three frames to this gene
			if(!_genes.get(geneID).containsTranscript(transcriptID+"_1")){
				_genes.get(geneID).addTranscript(transcriptID+"_1", chromosome, strand, transcriptBiotype, geneID, transcriptSymbol);
				_genes.get(geneID).addTranscript(transcriptID+"_2", chromosome, strand, transcriptBiotype, geneID, transcriptSymbol);
				_genes.get(geneID).addTranscript(transcriptID+"_3", chromosome, strand, transcriptBiotype, geneID, transcriptSymbol);

				// Set transcript biotype
				//String tmpAtt = trimAttribute(bits[19].trim());
				//genes.get(geneID).getTranscript(transcriptID+"_1").setBiotype(tmpAtt);
				//genes.get(geneID).getTranscript(transcriptID+"_2").setBiotype(tmpAtt);
				//genes.get(geneID).getTranscript(transcriptID+"_3").setBiotype(tmpAtt);

				/*
				if(_genes.get(geneID).getgeneSymbol().equals(_printVerboseForGeneName)){
					System.err.println(transcriptID+"_1: "+_isoformPeptideCounts.get(transcriptID+"_1"));
					System.err.println(transcriptID+"_2: "+_isoformPeptideCounts.get(transcriptID+"_2"));
					System.err.println(transcriptID+"_3: "+_isoformPeptideCounts.get(transcriptID+"_3"));
				}

				//System.out.println(transcriptID+"_1");
				if(_isoformPeptideCounts.containsKey(transcriptID+"_1"))
					tmp_peptideCount_1 = _isoformPeptideCounts.get(transcriptID+"_1");
				if(_isoformPeptideCounts.containsKey(transcriptID+"_2"))
					tmp_peptideCount_2 = _isoformPeptideCounts.get(transcriptID+"_2");
				if(_isoformPeptideCounts.containsKey(transcriptID+"_3"))
					tmp_peptideCount_3 = _isoformPeptideCounts.get(transcriptID+"_3");
				 */

				// Add gene-transcript map
				if(!transcriptID_2_geneID.containsKey(transcriptID+"_1")){
					transcriptID_2_geneID.put(transcriptID+"_1", geneID);
					transcriptID_2_geneID.put(transcriptID+"_2", geneID);
					transcriptID_2_geneID.put(transcriptID+"_3", geneID);
				}
			}

			// Add transcript or CDS length in nucleotides
			if(featureType.equals("CDS")){
				_genes.get(geneID).getTranscript(transcriptID+"_1").addCDS(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
				_genes.get(geneID).getTranscript(transcriptID+"_2").addCDS(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
				_genes.get(geneID).getTranscript(transcriptID+"_3").addCDS(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
				_genes.get(geneID).setProteinCoding(true);
			}else if(featureType.equals("exon")){
				_genes.get(geneID).getTranscript(transcriptID+"_1").addExon(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
				_genes.get(geneID).getTranscript(transcriptID+"_2").addExon(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
				_genes.get(geneID).getTranscript(transcriptID+"_3").addExon(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
			}


			/*
			// Reset exon lengths based on in-silico digest
			_genes.get(geneID).getTranscript(transcriptID+"_1").setExonLength(tmp_peptideCount_1);
			_genes.get(geneID).getTranscript(transcriptID+"_2").setExonLength(tmp_peptideCount_2);
			_genes.get(geneID).getTranscript(transcriptID+"_3").setExonLength(tmp_peptideCount_3);
			 //
			 // TODO: calculate peptides from the CDS sequence, rather than from the whole transcript
			 //
			_genes.get(geneID).getTranscript(transcriptID+"_1").setCodingExonLength(tmp_peptideCount_1);
			_genes.get(geneID).getTranscript(transcriptID+"_2").setCodingExonLength(tmp_peptideCount_2);
			_genes.get(geneID).getTranscript(transcriptID+"_3").setCodingExonLength(tmp_peptideCount_3);


			if(_genes.get(geneID).getgeneSymbol().equals(_printVerboseForGeneName)){
				_genes.get(geneID).getTranscript(transcriptID+"_1").printDebugInfo();
				_genes.get(geneID).getTranscript(transcriptID+"_2").printDebugInfo();
				_genes.get(geneID).getTranscript(transcriptID+"_3").printDebugInfo();
			}

			 */
			// Add feature length to the transcript/CDS length
			//int length = Integer.valueOf(bits[4]).intValue()-Integer.valueOf(bits[3]).intValue();
			//if(bits[2].equals("CDS"))
			//	genes.get(geneID).getTranscript(transcriptID).addCodingExonLength(length);
			//else
			//	genes.get(geneID).getTranscript(transcriptID).addExonLength(length);
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
		IO_utils.printLineErr("Reading expression data for priors: "+inputFile.getAbsolutePath());
		BufferedReader in = new BufferedReader(new FileReader(inputFile));
		String line = "";
		String[] bits;


		// determine from the file header if this is an RNA-seq (eXpress) output or a Thunder output form the footprint EM
		String inputDiscriminator = in.readLine().split("\t")[0];

		while((line=in.readLine()) != null){
			bits = line.split("\t");

			//transcriptQuants.put(bits[1].trim(), Double.valueOf(bits[14].trim()).doubleValue());
			String transcriptID = bits[1].trim();
			if(this.transcriptID_2_geneID.containsKey(transcriptID+"_1")){
				double tmpTPM = 0.0;
				if(inputDiscriminator.equals("bundle_id")){		// this is an eXpress file
					tmpTPM = Double.valueOf(bits[14].trim()).doubleValue();
				}else if(inputDiscriminator.equals("geneID")){		// this is a Thunder footprintEM file
					// "geneID\ttranscriptID\ttranscriptLength_exons\ttranscriptLength_CDS\tnFootprintsMappedToThisGene\tnReadsMappedToThisTranscript\thasFlatPrior\tprior\tflag\tnIterations\teffectiveReadCount\treadDensity_transcript\treadDensity_CDS\ttranscriptFractionAfterFootprints\ttranscriptBiotype");
					// use the transcript fraction as the TPM, this will get translated to the same prior
					tmpTPM = Double.valueOf(bits[13].trim()).doubleValue();
				}else{
					// Dunno!
				}
				this._genes.get(this.transcriptID_2_geneID.get(transcriptID+"_1")).getTranscript(transcriptID+"_1").setTPM(tmpTPM);
				this._genes.get(this.transcriptID_2_geneID.get(transcriptID+"_2")).getTranscript(transcriptID+"_2").setTPM(tmpTPM);
				this._genes.get(this.transcriptID_2_geneID.get(transcriptID+"_3")).getTranscript(transcriptID+"_3").setTPM(tmpTPM);

			}
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


	/**
	 * 
	 * @param maxIterations
	 * @param converganceDistance
	 */
	public void doEM(int maxIterations, double converganceDistance, boolean outputAll, boolean restrictToCDS) throws IOException{
		IO_utils.printLineErr("Running EM (maxIterations: "+maxIterations+"), writing to: "+this._outputPrefix+".exprs");
		Iterator<String> iterator = _genes.keySet().iterator();
		int totalGenes = _genes.size();
		int count = 0;
		int percent = 0;
		String geneID = "";
		EM_Core_optimised emCore;

		// variables to store the flags for summary at the end
		int count_genes_considered = 0, count_flag_converged = 0, count_flag_unconverged = 0, count_flag_error = 0; 

		PrintWriter out = new PrintWriter(new FileWriter(this._outputPrefix+".exprs"));

		if(_ms1Intensities.size() == 0)
			out.println("geneID\ttranscriptID\tframe\tnPeptidesMappedToThisGene\tnumberObservedPeptides\tnumberObservablePeptides\thasFlatPrior\tprior\tflag\tnIterations\teffectivePeptideCount\ttranscriptFractionAfterPeptides\tgeneSymbol\ttranscriptName\ttranscriptBiotype");
		else
			out.println("geneID\ttranscriptID\tframe\tnPeptidesMappedToThisGene\tnumberObservedPeptides\tnumberObservablePeptides\thasFlatPrior\tprior\tflag\tnIterations\taveragePrecursorIntensity\ttranscriptFractionAfterPeptides\tgeneSymbol\ttranscriptName\ttranscriptBiotype");

		while(iterator.hasNext()){
			geneID = iterator.next();

			if(_genes.get(geneID).getgeneSymbol().equals(_printVerboseForGeneName)){
				System.err.println("\n>doEM");
				_genes.get(geneID).printInfo();

				//System.err.println("\n>doEM");
				System.err.print("\n\nGeneID\tReadCountForThisGene\tTranscriptID\tObservablePeptides\tReadCountForThisTranscript\n");
			}

			// Change the effective exon length based on the number of mapped peptides (use # theoretical peptides when the in silico digest works well enough...)
			if(_genes.get(geneID).getReadCounts().size() > 0){
				for(int i=0;i<_genes.get(geneID).getTranscriptIDs().size();i++){
					String thisTranscriptID = _genes.get(geneID).getTranscriptIDs().get(i);
					if(_genes.get(geneID).getgeneSymbol().equals(_printVerboseForGeneName)){
						System.err.println(geneID+"\t"+_genes.get(geneID).getReadCounts().size()+"\t"+thisTranscriptID+"\t"+_genes.get(geneID).getTranscript(thisTranscriptID).getExonLength()+" / "+_genes.get(geneID).getTranscript(thisTranscriptID).getCodingExonLength()+"\t"+_genes.get(geneID).getTranscriptReadCounts(thisTranscriptID));
						//_genes.get(geneID).getTranscript(thisTranscriptID).printDebugInfo();
					}

					//_genes.get(geneID).getTranscript(thisTranscriptID).setExonLength(_genes.get(geneID).getReadCounts().size());
					//_genes.get(geneID).getTranscript(thisTranscriptID).setCodingExonLength(_genes.get(geneID).getReadCounts().size());
				}
			}

			// Run the EM
			emCore = new EM_Core_optimised(_genes.get(geneID), restrictToCDS);
			emCore.setReadLength(0);

			if(_genes.get(geneID).getgeneSymbol().equals(_printVerboseForGeneName)){
				emCore.setVerbose(true);
				//_genes.get(geneID).printInfo();
			}

			EM_Result res = emCore.runEM(maxIterations, converganceDistance);//, out);

			// Format and print EM results
			if(outputAll  ||  res.getGene().getReadCounts().size() > 0){

				// convergence status
				count_genes_considered ++;
				if(res.getFlag() == "1"){
					count_flag_converged ++;
				}else if(res.getFlag() == "0"){
					count_flag_unconverged ++;
				}else if(res.getFlag() == "-1"){
					count_flag_error ++;
				}

				for(int i=0;i<res.getGene().getTranscriptIDs().size();i++){
					String thisTranscriptID = res.getGene().getTranscriptIDs().get(i);
					int nPeptides = res.getGene().getTranscriptReadCounts(thisTranscriptID);
					double tmp_avgPrecursorIntensity = 0.0;
					if(nPeptides > 0){
						if(_ms1Intensities.size() == 0)
							tmp_avgPrecursorIntensity = res.getFinalEffectiveReadCount().get(thisTranscriptID).doubleValue()+0.0;
						else
							tmp_avgPrecursorIntensity = (res.getFinalEffectiveReadCount().get(thisTranscriptID).doubleValue()+0.0) / nPeptides;
					}
					int tmp_coverage = res.getGene().getTranscript(thisTranscriptID).getCodingExonLength();
					double tmp_finalLikelihood = res.getFinalTranscriptLikelihoods().get(thisTranscriptID).doubleValue();
					//out.printf(geneID+"\t"+thisTranscriptID+"\t"+res.getGene().hasFlatPrior()+"\t%e\t"+res.getFlag()+"\t"+res.getNIterations()+"\t%f"+nPeptides+"\t%e\t%e\n", res.getGene().getPriors().get(thisTranscriptID), tmp_avgPrecursorIntensity, tmp_coverage, tmp_finalLikelihood);
					String[] transcriptIDbits = thisTranscriptID.split("_");
					out.printf(geneID+"\t"+
							transcriptIDbits[0]+"\t"+
							transcriptIDbits[1]+"\t"+
							res.getGene().getReadCounts().size()+"\t"+
							nPeptides+"\t"+
							tmp_coverage+"\t"+
							res.getGene().hasFlatPrior()+"\t%e\t"+
							res.getFlag()+"\t"+
							res.getNIterations()+"\t%f\t%e\t"+
							res.getGene().getgeneSymbol()+"\t"+
							res.getGene().getTranscript(thisTranscriptID).getTranscriptSymbol()+"\t"+
							res.getGene().getTranscript(thisTranscriptID).getTranscriptBiotype()+"\n",
							res.getGene().getPriors().get(thisTranscriptID), tmp_avgPrecursorIntensity, tmp_finalLikelihood);
				}
			}

			// Update the status bar
			count++;
			if(Math.round((count*100.0)/totalGenes) > percent){
				percent = (int)Math.round((count*100.0)/totalGenes);
				IO_utils.printProgressBar(percent);
			}
		}
		System.err.println("");
		out.flush();
		out.close();

		//IO_utils.printLineErr("  "+count_genes_considered+" genes with valid peptide(s) ("+Math.round(count_genes_considered*10000/count_genes_considered)/100.0+"% of total annotations)");
		IO_utils.printLineErr("  "+count_flag_converged+" genes converged successfully ("+Math.round(count_flag_converged*10000/count_genes_considered)/100.0+"%)");
		IO_utils.printLineErr("  "+count_flag_unconverged+" genes failed to converge ("+Math.round(count_flag_unconverged*10000/count_genes_considered)/100.0+"%)");
		IO_utils.printLineErr("  "+count_flag_error+" genes resulted in an error ("+Math.round(count_flag_error*10000/count_genes_considered)/100.0+"%)");

		IO_utils.printLineErr("Done.");
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withLongOpt("annotation").withArgName("path").hasArg().withDescription("Path to the GTF file containing the transcript/gene relationship").create(MIBAT.OPT_PATH_ANNOTATION));
		options.addOption(OptionBuilder.withLongOpt("fasta").withArgName("path").hasArg().withDescription("Path to the fasta file containing [AA] isoform sequences so that we can predict expected peptide products").create("s"));
		options.addOption(OptionBuilder.withLongOpt("priors").withArgName("path").hasArg().withDescription("[optional] Text file containing RNA-seq transcript expression quantifications from eXpress").create("e"));
		options.addOption(OptionBuilder.withLongOpt("spectra_mapped").withArgName("path").hasArg().withDescription("X!Tandem output XML file containing spectra alignments.  Multiple input files can be specified using a comma separator and no space, these will be merged and analysed together.").create("f"));
		options.addOption(OptionBuilder.withArgName("peptidesPath").hasArg().withDescription("Max Quant output table containing peptide alignments.").create("maxquant"));
		options.addOption(OptionBuilder.withLongOpt("sampleID").withArgName("string").hasArg().withDescription("ID of the sample to extract from the maxquant peptides table").create("id"));
		options.addOption(OptionBuilder.withLongOpt("spectra_raw").withArgName("path").hasArg().withDescription("[optional] mzXML file containing containing raw spectra peaks (required for quantification)").create("mzxml"));		
		options.addOption(OptionBuilder.withLongOpt("outputPrefix").withArgName("outputPath").hasArg().withDescription("Path to which to output results").create("o"));
		options.addOption(OptionBuilder.withDescription("Count only peptides in the CDS (where available)").create("cds"));
		options.addOption(OptionBuilder.withArgName("int").hasArg().withDescription("Maximum number of EM iterations [default: 10]").create("N"));
		options.addOption(OptionBuilder.withArgName("double").hasArg().withDescription("Criteria for EM convergence [default: 0.0001]").create("c"));
		options.addOption(OptionBuilder.withArgName("double").hasArg().withDescription("[default: 0.05] Based on hits to the reverse decoy DB, define an acceptable false discovery rate (FDR)").create("fdr"));
		options.addOption(OptionBuilder.withDescription("Write all genes/isoforms, even those with no observed peptide spectra").create("v"));
		return options;
	}




	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {

		/*
		args = new String[]{"IsoformEM_Proteomics", 
				//"-f", "/Users/robk/Desktop/EM_TEST/ProteoData/allSpectra_C8.txt",
				"-f", "/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Proteomics/HEK.C8_CEX.Orbi.gencode21.TandemOutput.xml",
				//"-f", "/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Proteomics/HEK.C8.Orbi.gencode21.TandemOutput.xml",
				//"-mzxml", "/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Proteomics/EOT12-0395.HEK.C8.Orbi.mzxml",
				"-a", "/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v21.annotation_noSelenocysteine.gtf",
				//"-a", "/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Proteomics/ProteomicEM/tmp.ann",
				//"-e", "/Users/robk/Box Sync/Work/HEK293_RNAseq/STAR_eXpress/A1-Total_2x75/results.xprs",
				//"-e", "/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/A1_totalRNA/results.xprs",
				"-e", "/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Footprinting/FootprintEM/{sample=FP_B3}{prior=A2_L10aRNA}{filter=cds}{maxIt=1000}{frame=use}.exprs",
				"-s", "/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v21.annotation.protein.fa",
				"--outputPrefix","/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Proteomics/ProteomicEM/TEST_TEST_TEST",
				"--cds"};

		//"-o", "/Users/robk/Desktop/EM_TEST/TEST_NEW",
		//"-DEV"};
		 */

		/*
		args = new String[]{"IsoformEM_Proteomics",
				//"--maxquant","/Users/robk/Downloads/MQ_TEST_ALL.tsv",
				"--maxquant","/Users/robk/Dropbox/MaxQuant output CEGs/peptides.txt",
				"--sampleID","HSB105_AMY",
				"-a", "/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v21.annotation_noSelenocysteine.gtf",
				"-e", "/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Footprinting/FootprintEM/{sample=FP_B3}{prior=A2_L10aRNA}{filter=cds}{maxIt=1000}{frame=use}.exprs",
				"-s", "/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v21.annotation.protein.fa",
				"--outputPrefix","/Users/robk/Downloads/TEST_TEST_TEST"};
		 */
		
		CommandLine cmdArgs = MIBAT.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption(MIBAT.OPT_PATH_ANNOTATION)  &&  (cmdArgs.hasOption("f") || cmdArgs.hasOption("maxquant"))  &&  cmdArgs.hasOption("s") && cmdArgs.hasOption("o")){
			System.err.println();

			ProteomicEM engine = new ProteomicEM(cmdArgs.getOptionValue("o"));



			// for testing!
			//engine._printVerboseForGeneName = "COG8";
			//engine._printVerboseForGeneName = "SARS";

			// Clean up
			//if(cmdArgs.hasOption("v")){
			//engine.removeFile(new File(tmp_path));
			//	engine.verbose = true;
			//}

			//String output_path = cmdArgs.getOptionValue("o");
			//File output_priors = new File(output_path+"/priors.txt");
			//File output_alignments = new File(tmp_path+"/alignments.txt");


			// Create the output directory and temp dir if they don't already exist
			//engine.makeOutputDirectories(output_path);

			// read the number of expected peptides per isoform per frame
			engine.readDigestInfo(new File(cmdArgs.getOptionValue("s")));

			// Read the transcripts in the GTF
			engine.readGTF(new File(cmdArgs.getOptionValue(MIBAT.OPT_PATH_ANNOTATION)));

			// Read the eXpress transcript quantifications or choose flat priors
			if(cmdArgs.hasOption("e")){
				engine.readTranscriptExpressions(new File(cmdArgs.getOptionValue("e")));
			}

			// Convert the transcript expressions to priors and output
			engine.initialisePriors();

			// Set desired FDR for spectra
			double fdrMax = 0.05;
			if(cmdArgs.hasOption("fdr"))
				fdrMax = Double.valueOf(cmdArgs.getOptionValue("fdr")).doubleValue();

			// if this is a maxquant file read it, else expect Tandem XML
			if(cmdArgs.hasOption("maxquant")){
				//engine.readPeptides_maxquant(cmdArgs.getOptionValue("maxquant"), "");
				engine.readPeptides_maxquant(cmdArgs.getOptionValue("maxquant"), cmdArgs.getOptionValue("sampleID"));
			}else{
				// If provided an mzXML file containing precursor intensities, read it now
				if(cmdArgs.hasOption("mzxml"))
					engine.readMZXML(new File(cmdArgs.getOptionValue("mzxml")));

				// Parse the MS/MS alignments
				engine.readSpectra_tandem(cmdArgs.getOptionValue("f"), fdrMax, cmdArgs.hasOption("cds"));
			}

			// set the number of iterations to cap the EM
			int maxIterations = 10;
			if(cmdArgs.hasOption("N"))
				maxIterations = Integer.valueOf(cmdArgs.getOptionValue("N")).intValue();

			// set the convergence criteria for the EM
			double convergenceDistance = 1.0/10000.0;
			if(cmdArgs.hasOption("c"))
				convergenceDistance = Double.valueOf(cmdArgs.getOptionValue("c")).doubleValue();

			boolean outputAll = false;
			if(cmdArgs.hasOption("v"))
				outputAll = true;

			// set all alignment weights to 1 (as there is nothing yet like the footprint frames that can modify this for MSMS)
			//engine.resetAlignmentWeights();

			// Run the EM!
			engine.doEM(maxIterations, convergenceDistance, outputAll, true);

			IO_utils.printLineErr("All Done!");

		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.setWidth(200);
			formatter.printHelp(MIBAT.THUNDER_EXE_COMMAND+" IsoformEM_Proteomics", getCmdLineOptions());
			System.err.println();
		}


	}

	/*public String getTime(){
		return((new SimpleDateFormat("yyyy/MM/dd HH:mm:ss")).format(Calendar.getInstance().getTime()));
	}*/




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



