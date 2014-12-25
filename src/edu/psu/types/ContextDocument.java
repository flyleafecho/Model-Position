package edu.psu.types;

import java.util.*;
import gnu.trove.*;
import edu.psu.util.RemoveStopWord;
import java.io.*;

/*
 * Data structure
 * Initilization
 * Sampling
 * Results
 */

public class ContextDocument extends Document {
	
	// in its constructor
	public TIntIntHashMap citationSet;
	//add two variables for different parts
	public TIntIntHashMap[] citationSetByPart;
	
	public TIntIntHashMap citationTopicAssignment;
	//add citation topic assignment for different part
	public int[][] citationTAByPart;

	public ContextDocument() {
		partNums = 2;
		citationSet = new TIntIntHashMap();
		citationSetByPart = new TIntIntHashMap[partNums];
	}
	
	public void add_milestone(String name, int docid, String context, ArrayList<String> Docs,ArrayList<String> citations)
	{
		int label = 0;
		Vector<String> contentVector = new Vector<String>();
		Vector<Citation> citationVector = new Vector<Citation>();
		ArrayList<String> docContent=new ArrayList<String>();
		int DocIndex=Docs.indexOf(name);
		//Format of Str : #cid:#part\t#cid:#part.....
		String[] citation_part = citations.get(DocIndex).split("\\t");	
		for(int i=0;i<citation_part.length;i++)
		{
			String citation = citation_part[i].split(":")[0];
			int part = Integer.parseInt(citation_part[i].split(":")[1]);
			citationVector.add(new Citation(citation,part)) ;
			Corpus.citationAlphabet.lookupIndex(citation);
			citationSet.adjustOrPutValue(Corpus.citationAlphabet.lookupIndex(citation), 1,1);
			citationSetByPart[part].adjustOrPutValue(Corpus.citationAlphabet.lookupIndex(citation), 1, 1);			
		}
			
		String [] words=context.split("\\W+");
		for(String word : words){
			String tmp = word.trim();
			if(tmp.trim().length()>2){
				docContent.add(tmp);
			}
		}	
		
		RemoveStopWord removeStopWords =new RemoveStopWord(new File("D:/WorkPaceEclipse_2/citation_LDA-master/src/stoplists/en.txt"), "UTF-8",docContent.toArray(new String[docContent.size()]), false, false);
	    String[] removeStopTokens=removeStopWords.RemoveStopWords();
	    for(int k=0;k<removeStopTokens.length;k++)
	    {
	    	if(!removeStopTokens[k].matches(".*\\d+.* | \\W+ |\\d+"))
	    	{
	    		if(!removeStopTokens[k].replace("[\\p{Punct}\\p{Space}]", "").trim().equals(""))
	    		{
	    			contentVector.add(removeStopTokens[k].replace("[\\p{Punct}\\p{Space}]", ""));
	    		}
	    	}	

	    }
		
		ProcessCSV(name, label, docid, contentVector, citationVector);	
	}

}