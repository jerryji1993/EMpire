package expectationMaximisation;

public class EM_Alignment{

	private String _readID, _transcriptID;
	//private double _alignmentWeight = 1.0/3.0;
	private double _alignmentWeight = 0.0;
	
	public EM_Alignment(String readID, String transcriptID){
		_readID = readID;
		_transcriptID = transcriptID;
	}
	
	public EM_Alignment(String readID, String transcriptID, double weight){
		_readID = readID;
		_transcriptID = transcriptID;
		_alignmentWeight = weight;
	}
	
	public String getReadID(){ return _readID; }
	public String getTranscriptID(){ return _transcriptID; }
	public double getAlignmentWeight(){ return _alignmentWeight; }
	
	public void setAlignmentWeight(double alignmentWeight){ _alignmentWeight = alignmentWeight; }
	
}
