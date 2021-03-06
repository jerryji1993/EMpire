package expectationMaximisation;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import objects.Transcript;

public class EM_Core_optimised {

	private boolean _verbose = false;
	public void setVerbose(boolean verbose){ _verbose = verbose; }
	
	private EM_Gene _gene;
	private ArrayList<String> _transcriptIDs;
	//private boolean _hasFlatPrior = false;

	double[] _priors;//, _finalTranscriptRatios, _finalAbundances;
	int[] _transcriptLengths;

	//private HashMap<String, ArrayList<String>> _reads = new HashMap<String, ArrayList<String>>();
	private Object[] _readIDs;

	private static final double _TINY = 1.0/10000.0;



	public EM_Core_optimised(EM_Gene gene, boolean useCDS){
		_gene = gene;
		_transcriptIDs = _gene.getTranscriptIDs();
		Collections.sort(_transcriptIDs);

		Transcript tmp;

		_priors = new double[_transcriptIDs.size()];
		_transcriptLengths = new int[_transcriptIDs.size()];

		for(int i=0;i<_transcriptIDs.size();i++){
			tmp = _gene.getTranscript(_transcriptIDs.get(i));
			//System.out.println("transcript: "+tmp.getID());
			_priors[i] = tmp.getPrior();
			if(useCDS)
				_transcriptLengths[i] = tmp.getTotalCodingExonLength();
			else
				_transcriptLengths[i] = tmp.getTotalExonLength();

			//if(!_hasFlatPrior)
			//_hasFlatPrior = !tmp.hasRNAseq();
		}
		//_reads = _gene.getReads();
		//_readIDs = _gene.keySet().toArray();
		_readIDs = _gene.getReadIDs().toArray();
	}

	/**
	 * 
	 * @param keepRead2isoform
	 */
//	private void filterRead2TranscriptAlignments(HashMap<String, HashSet<String>> keepRead2isoform){
//
//		Iterator<String> reads = _reads.keySet().iterator();
//		// make a list of reads to remove that are not in frame with any isoform
//		HashSet<String> readsToRemove = new HashSet<String>();
//
//		while(reads.hasNext()){
//			String thisReadID = reads.next();
//
//			// make a list of isoforms to remove that are not in frame with this read
//			HashSet<String> isoformsToRemove = new HashSet<String>();
//
//			// if this read does not appear in the frame list, remove it
//			if(!keepRead2isoform.containsKey(thisReadID)){
//				readsToRemove.add(thisReadID);
//			}
//			else{
//				// add isoforms to remove
//				Iterator<String> isoforms = _reads.get(thisReadID).iterator();
//				while(isoforms.hasNext()){
//					String thisIsoform = isoforms.next();
//
//					if(!keepRead2isoform.get(thisReadID).contains(thisIsoform))
//						isoformsToRemove.add(thisIsoform);
//				}
//
//				// remove the isoforms that are out of frame
//				for(String removeIsoform: isoformsToRemove){
//					_reads.get(thisReadID).remove(removeIsoform);
//					_gene.removeIsoformReadCount(removeIsoform);
//				}
//
//				// if there are no in frame alignments for this read, add it to the list to be removed
//				if(_reads.get(thisReadID).size() == 0)
//					readsToRemove.add(thisReadID);
//			}
//		}
//
//		// remove the reads that are out of frame with all isoforms
//		for(String removeRead: readsToRemove){
//			_reads.remove(removeRead);
//			_gene.removeRead(removeRead);
//		}
//
//		// re-compute array of reads
//		_readIDs = _reads.keySet().toArray();
//	}


	public EM_Result runEM(int maxIterations, double converganceDistance){//, PrintWriter out_exprs){

		
		if(_verbose){
			System.err.println("\n>runEM");
			System.err.println("ReadID\ttranscriptID\tweight");
			for(String read: _gene.getReadIDs()){
				for(String transcript: _gene.getTranscriptsAlignedToRead(read))
					System.err.println(read+"\t"+transcript+"\t"+_gene.getAlignment(read, transcript).getAlignmentWeight());
			}
		}
		
		HashMap<String,BigDecimal> tmpReadCount = new HashMap<String,BigDecimal>();
		HashMap<String,BigDecimal> tmpPriors = new HashMap<String,BigDecimal>();
		HashMap<String,BigDecimal> tmpFinalLikelihoods = new HashMap<String,BigDecimal>();

		if(_gene.transcriptCount() == 1){									// single isoform gene!
			//out_exprs.printf(_gene.getID()+"\t"+_transcriptIDs.get(0)+"\t"+_hasFlatPrior+"\t%e\t1\t0\t%f\t%e\t%e\n", _priors[0]+0.0, _reads.size()+0.0, _reads.size()/(_transcriptLengths[0]-_readLength+1.0), _priors[0]+0.0);
			tmpReadCount.put(_transcriptIDs.get(0), new BigDecimal(_gene.getReadCount()+0.0));
			tmpPriors.put(_transcriptIDs.get(0), new BigDecimal(_priors[0]+0.0));
			return(new EM_Result(_gene, 0, "1", tmpReadCount, tmpPriors));
		}else if(_gene.getReadCount() == 0){										// no reads!
			for(int k=0;k<_transcriptIDs.size();k++){
				//out_exprs.printf(_gene.getID()+"\t"+_transcriptIDs.get(k)+"\t"+_hasFlatPrior+"\t%e\t1\tNA\t%f\t%e\t%e\n", _priors[k]+0.0, 0.0, 0.0, _priors[k]+0.0);
				tmpReadCount.put(_transcriptIDs.get(k), new BigDecimal(0.0));
				tmpPriors.put(_transcriptIDs.get(k), new BigDecimal(_priors[k]+0.0));
			}
			return(new EM_Result(_gene, 0, "NA", tmpReadCount, tmpPriors));
		}else{																// multi-isoform gene with read support, do EM!

			_mu = _priors;
			double[] mu_new;
			double dist = 0.0;
			String flag = "0";
			
			if(_verbose){
				System.err.print(_gene.getID()+"\tIDs=");
				for(int k=0;k<_transcriptIDs.size();k++)
					System.err.print(_transcriptIDs.get(k)+",");
				System.err.print("\n"+_gene.getID()+"\tpriors=");
				for(int k=0;k<_priors.length;k++)
					System.err.print(_priors[k]+",");
				System.err.println();
				
			}
			
			int count = 0;
			for(int i=0;i<maxIterations;i++){

				// E-step
				expectationStep();

				// M-step
				mu_new = maximisationStep();

				// Distance calculation
				dist = this.calcConvergence(_mu, mu_new);

				if(_verbose){
					//System.err.print(_gene.getID()+"\tk="+i+"\tmu=");
					//for(int k=0;k<_mu.length;k++)
					//	System.err.print(_mu[k]+",");
					System.err.print(_gene.getID()+"\tk="+i+"\tmu_new=");
					for(int k=0;k<_mu.length;k++)
						System.err.print(mu_new[k]+",");
					System.err.println("\tdis="+dist);
				}

				// Update
				_mu = mu_new;
				count ++;

				// Stop?
				if(dist < converganceDistance){
					flag = "1";
					break;
				}
				double sum = 0.0;
				for(int k=0;k<mu_new.length;k++)
					sum += mu_new[k];
				if(sum <= 0.0){
					flag = "-1";
					break;
				}


			}

			// Finished with this gene
			for(int k=0;k<_transcriptIDs.size();k++){
				//out_exprs.printf(_gene.getID()+"\t"+_transcriptIDs.get(k)+"\t"+_hasFlatPrior+"\t%e\t"+flag+"\t"+count+"\t%f\t%e\t%e\n", _priors[k]+0.0, _readCounts[k]+0.0, _readCounts[k]/(_transcriptLengths[k]-_readLength+1.0), _mu[k]+0.0);
				tmpReadCount.put(_transcriptIDs.get(k), new BigDecimal(_readCounts[k]+0.0));
				tmpPriors.put(_transcriptIDs.get(k), new BigDecimal(_priors[k]+0.0));
				tmpFinalLikelihoods.put(_transcriptIDs.get(k), new BigDecimal(_mu[k]+0.0));
			}
			return(new EM_Result(_gene, count, flag, tmpReadCount, tmpFinalLikelihoods));
		}
	}



	public HashMap<String, String> sampleReadAssignments(HashMap<String, String> readMap){
		Iterator<String> reads = _gene.getReadIDs().iterator();
		String thisRead;
		//ArrayList<String> theseTranscripts;
		//HashMap<String, String> readMap = new HashMap<String, String>();

		while(reads.hasNext()){
			thisRead = reads.next();
			//theseTranscripts = _reads.get(thisRead);
			//double[] transcriptFracs = new double[_reads.get(thisRead).size()];
			double[] transcriptFracs = new double[_transcriptIDs.size()];
			double transcriptSum = 0.0;

			if(_transcriptIDs.size() == 1){
				transcriptFracs[0] = 1.0;
				transcriptSum = 1.0;
			}else{
				for(int k=0;k<_transcriptIDs.size();k++){
					String thisTranscriptID = _transcriptIDs.get(k);
					//if(theseTranscripts.contains(thisTranscriptID))
					if(_gene.hasAlignment(thisRead, thisTranscriptID)){
						transcriptFracs[k] = _readCounts[k];
						transcriptSum += _readCounts[k]+0.0;
					}
					else{
						transcriptFracs[k] = 0.0;
					}
				}
			}

			for(int k=0;k<_transcriptIDs.size();k++){
				transcriptFracs[k] = transcriptFracs[k] / (transcriptSum+0.0);
			}

			double rand = Math.random();
			double cumsum = 0.0;
			for(int k=0;k<_transcriptIDs.size();k++){
				cumsum += transcriptFracs[k];
				if(cumsum >= rand){
					readMap.put(thisRead, _transcriptIDs.get(k));
					break;
				}
			}
		}
		return(readMap);
	}


	private int _readLength = 29;
	public void setReadLength(int readLength){ this._readLength = readLength; }

	private double[] mu2tau(double[] priors){
		double[] tau = new double[priors.length];

		if(priors.length != _transcriptLengths.length){
			System.err.println("mu2tau: input vectors are not of the same length\n");
			System.exit(1);
		}

		double tauSum = 0.0;
		for(int k=0;k<priors.length;k++){
			tau[k] = priors[k] * Math.max(0, _transcriptLengths[k]-_readLength);
			tauSum += tau[k]; 
		}

		if(Math.abs(tauSum - _TINY) < 0){
			System.err.println("mu_2_tau: sum of prior is too small, please check parameters\n");
			System.exit(1);
		}

		for(int k=0;k<priors.length;k++){
			tau[k] = tau[k]/tauSum;
		}

		return(tau);
	}


	private double[] allocateRead(String readID, double[] tau, double[] readCounts){
		double[] addVal = new double[tau.length];
		double sumVal = 0.0;
		for(int k=0;k<_transcriptIDs.size();k++){
			String thisTranscriptID = _transcriptIDs.get(k);
			if(_gene.hasAlignment(readID, thisTranscriptID)){
				//_reads.get(readID).get(thisTranscriptID);
				//addVal[k] = tau[k];
				
				addVal[k] = tau[k] * _gene.getAlignment(readID, thisTranscriptID).getAlignmentWeight();
			}else{
				addVal[k] = 0.0;
			}
			sumVal += addVal[k];
		}

		if(sumVal > 0.0){
			//*^^^ split the read proportionally to the transcripts^^^*/
			for(int k=0;k<_transcriptIDs.size();k++){
				//readCounts[k] = readCounts[k] + addVal[k]/sumVal;
				readCounts[k] = readCounts[k] + ((addVal[k]/sumVal) * _gene.getReadCount(readID));
			}
		}

		return(readCounts);
	}


	// mu are the priors
	private double[] _mu, _tau, _readCounts;

	/**
	 * 
	 */
	public void expectationStep(){

		_tau = mu2tau(_mu);
		_readCounts = new double[_tau.length];

		for(int k=0; k<_readIDs.length; k++){
			_readCounts = allocateRead((String)_readIDs[k], _tau, _readCounts);
		}
	}

	/**
	 * 
	 * @return
	 */
	public double[] maximisationStep(){
		double[] mu_new = new double[_mu.length];
		double s = 0.0;

		for(int k=0;k<_readCounts.length;k++){
			mu_new[k] = _readCounts[k]/Math.max(1, _transcriptLengths[k]-_readLength);
			s += mu_new[k];
		}
		if(s > 0.0){
			for(int k=0;k<_readCounts.length;k++){
				mu_new[k] = mu_new[k]/s;
			}
		}
		return(mu_new);
	}

	/**
	 * 
	 * @param mu_old
	 * @param mu_new
	 * @return
	 */
	public double calcConvergence(double[] mu_old, double[] mu_new){
		if(mu_old.length != mu_new.length){
			System.err.println("cover_distance: two transcripts expression percent not of the same length\n");
			System.exit(1);
		}
		double dis = 0.0;
		for(int k=0;k<mu_old.length;k++){
			dis += Math.abs(mu_old[k]-mu_new[k]);
		}
		return(dis);
	}


}
