package expectationMaximisation;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import objects.Transcript;

public class EM_Gene{
	private String _id, _geneSymbol;
	//private HashMap<String, EM_Transcript> _transcripts = new HashMap<String, EM_Transcript>();	
	private HashMap<String, Transcript> _transcripts = new HashMap<String, Transcript>();

	public EM_Gene(String id, String geneSymbol){
		_id = id;
		_geneSymbol = geneSymbol;
	}
	public String getID(){ return _id; }
	public String getgeneSymbol(){ return _geneSymbol; }

	public boolean containsTranscript(String id){ return(_transcripts.containsKey(id)); }
	//public void addTranscript(String id, String transcriptSymbol){ 
	//	_transcripts.put(id, new EM_Transcript(id, transcriptSymbol));
	//}
	
	public void addTranscript(String id, String chromosome, String strand, String transcriptBiotype, String geneID, String transcriptSymbol){ 
		_transcripts.put(id, new Transcript(id, chromosome, strand, transcriptBiotype, geneID, transcriptSymbol));
	}


	//private boolean XX_restrictToCDS = false;
	//public void setRestrictToCDS(boolean val){ _restrictToCDS = val; }
	//public boolean IsRestrictedToCDS(){ return _restrictToCDS; }

	private boolean _hasCDS = false;
	public void setProteinCoding(boolean val){ _hasCDS = val; }
	public boolean hasCDS(){ return _hasCDS; }

	public int transcriptCount(){ return(_transcripts.size()); }
	//public EM_Transcript getTranscript(String id){ 
	//	return(_transcripts.get(id)); 
	//}
	public Transcript getTranscript(String id){ 
		return(_transcripts.get(id)); 
	}

	/**
	 * 
	 * @return
	 */
	public ArrayList<String> getTranscriptIDs(){
		ArrayList<String> ids = new ArrayList<String>();
		Iterator<String> it = _transcripts.keySet().iterator();
		while(it.hasNext()){
			ids.add(it.next());
		}
		Collections.sort(ids);
		return(ids);
	}


	/**
	 * 
	 */
	public void setPriors(){
		double totalTPM = 0.0;
		Iterator<String> it = _transcripts.keySet().iterator();
		while(it.hasNext()){
			totalTPM += _transcripts.get(it.next()).getTPM();
		}
		it = _transcripts.keySet().iterator();
		String id = "";
		while(it.hasNext()){
			id = it.next();
			if(totalTPM == 0.0){
				_transcripts.get(id).setPrior(1.0 / (_transcripts.size()+0.0));
				_hasFlatPrior = true;
			}else{
				_transcripts.get(id).setPrior(_transcripts.get(id).getTPM() / totalTPM);
				_hasFlatPrior = false;
			}
		}
	}

	private boolean _hasFlatPrior = true;
	public boolean hasFlatPrior(){ return _hasFlatPrior; }



	/**
	 * 
	 * @return
	 */
	public HashMap<String, BigDecimal> getPriors(){
		HashMap<String, BigDecimal> priors = new HashMap<String, BigDecimal>();
		Iterator<String> it = _transcripts.keySet().iterator();
		String id;
		while(it.hasNext()){
			id = it.next();
			priors.put(id, new BigDecimal(_transcripts.get(id).getPrior()));
		}
		return priors;
	}


	/**
	 * 
	 * @return
	 */
	public HashMap<String, Double> getAbundances(){
		HashMap<String, Double> expressions = new HashMap<String, Double>();
		Iterator<String> it = _transcripts.keySet().iterator();
		String id;
		while(it.hasNext()){
			id = it.next();
			expressions.put(id, _transcripts.get(id).getTPM());
		}
		return expressions;
	}

	/**
	 * 
	 * @return
	 */
	/*public HashMap<String, Integer> gXetLengths(){
		HashMap<String, Integer> lengths = new HashMap<String, Integer>();
		Iterator<String> it = _transcripts.keySet().iterator();
		String id;
		while(it.hasNext()){
			id = it.next();
			if(!_restrictToCDS)
				lengths.put(id, _transcripts.get(id).getExonLength());
			else
				lengths.put(id, _transcripts.get(id).getCodingExonLength());
		}
		return lengths;
	}*/


	// Stores read->transcript map:  HashMap<readID, ArrayList<transcriptIDs>>
	//private HashMap<String, ArrayList<String>> _reads = new HashMap<String, ArrayList<String>>();

	// HashMap<readID, HashMap<transcriptIDs, frameWeight>>	
	//private HashMap<String, HashMap<String, Double>> _reads2transcripts2weight = new HashMap<String, HashMap<String, Double>>();
	
	private HashMap<String, Long> _readCounts = new HashMap<String, Long>();
	//private HashMap<String, Integer> _transcriptReadCounts = new HashMap<String, Integer>();

	private HashMap<Integer, EM_Alignment> _index2alignment = new HashMap<Integer, EM_Alignment>();
	private HashMap<String, HashMap<String, Integer>> _reads2transcripts2alignmentIndex = new HashMap<String, HashMap<String, Integer>>();
	private HashMap<String, HashMap<String, Integer>> _transcripts2reads2alignmentIndex = new HashMap<String, HashMap<String, Integer>>();
	
	public boolean hasAlignment(String readID, String transcriptID){
		boolean result = false;
		if(_reads2transcripts2alignmentIndex.containsKey(readID))
			if(_reads2transcripts2alignmentIndex.get(readID).containsKey(transcriptID))
				result = true;
		return result;
	}
	
	public EM_Alignment getAlignment(String readID, String transcriptID){
		EM_Alignment result = null;
		if(hasAlignment(readID, transcriptID))
			result = _index2alignment.get(_reads2transcripts2alignmentIndex.get(readID).get(transcriptID));
		return result;
	}
	public double getAlignmentWeight(String readID, String transcriptID){
		return getAlignment(readID, transcriptID).getAlignmentWeight();
	}
	
	/**
	 * 
	 */
	public Set<String> getReadIDs(){ return _readCounts.keySet(); }
	

	
	public Set<String> getTranscriptsAlignedToRead(String readID){
		return _reads2transcripts2alignmentIndex.get(readID).keySet();
	}
	
	public int getReadCount(){ return _readCounts.size(); }
	public int getTranscriptReadCounts(String transcriptID){
		if(_transcripts2reads2alignmentIndex.containsKey(transcriptID))
			return _transcripts2reads2alignmentIndex.get(transcriptID).size();
		else
			return 0;
	}



	/**
	 * Add a read and alignments (meant for situations where a 'read' is more like a peptide in that it has some intensity value, rather than being an atomic unit
	 * @param readID
	 * @param hitsTranscripts
	 * @param readCount
	 */
	public void addRead(String readID, ArrayList<String> hitsTranscripts, long readCount){
		
		// if this read ID is not in the alignment index list
		if(!_reads2transcripts2alignmentIndex.containsKey(readID)){
			_reads2transcripts2alignmentIndex.put(readID, new HashMap<String, Integer>());
			_readCounts.put(readID, new Long(0));
		}
		
		// add read count
		_readCounts.put(readID, _readCounts.get(readID)+readCount);
		
		// loop through all transcripts to which this read maps
		for(String transcript: hitsTranscripts){
			
			// add this new alignments
			int thisAlignmentIndex = _index2alignment.size() + 1;
			_index2alignment.put(thisAlignmentIndex, new EM_Alignment(readID, transcript));
			
			// add to read -> transcript map
			_reads2transcripts2alignmentIndex.get(readID).put(transcript, thisAlignmentIndex);
			
			// add to transcript -> read map
			if(!_transcripts2reads2alignmentIndex.containsKey(transcript))
				_transcripts2reads2alignmentIndex.put(transcript, new HashMap<String, Integer>());
			_transcripts2reads2alignmentIndex.get(transcript).put(readID, thisAlignmentIndex);
		}
		
	}
	public void addRead(String readID, ArrayList<String> hitsTranscripts){
		addRead(readID, hitsTranscripts, 1);
	}


	/**
	 * If we're considering the footprint frame, add the weight 
	 * @param readID
	 * @param transcriptID
	 * @param weight
	 */
	public void addWeightToReadAlignment(String readID, String transcriptID, double weight){
		getAlignment(readID, transcriptID).setAlignmentWeight(weight);
	}
	
	
	/**
	 * Resets the alignment weights for all reads/transcripts in this gene (helpful if there are no weights that were added)
	 * @param newWeight
	 */
	public void resetAlignmentWeights(double newWeight){
		for(String transcriptID: _transcripts2reads2alignmentIndex.keySet()){
			for(String readID: _transcripts2reads2alignmentIndex.get(transcriptID).keySet())
				_index2alignment.get(_transcripts2reads2alignmentIndex.get(transcriptID).get(readID)).setAlignmentWeight(newWeight);
		}
	}
	
	
	//public HashMap<String, ArrayList<String>> getReads(){ return _reads; }
	public HashMap<String, Long> getReadCounts(){ return _readCounts; }
	public Long getReadCount(String readID){ return _readCounts.get(readID); }
	
	
	public void printInfo(){
		System.err.println(_id+"\t"+_geneSymbol);
		
		System.err.println("_readCounts:");
		for(String readID: _readCounts.keySet())
			System.err.println("\t> "+readID+" : "+_readCounts.get(readID));
		
		System.err.println("_transcripts2reads2alignmentIndex:");
		for(String transcriptID: _transcripts2reads2alignmentIndex.keySet()){
			for(String readID: _transcripts2reads2alignmentIndex.get(transcriptID).keySet()){
				System.err.println("\t> "+transcriptID+"\t"+readID+"\t"+_index2alignment.get(_transcripts2reads2alignmentIndex.get(transcriptID).get(readID)).getAlignmentWeight());
			}
		}
	}
}




