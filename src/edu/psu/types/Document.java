package edu.psu.types;

import java.util.*;

import gnu.trove.*;
import edu.psu.topic.Model;
import edu.psu.types.*;

import java.io.*;

public class Document implements Serializable {

	public LabelSequence topicAssignments = null; // sequence representation of topic assignments
	
	public int [] wordTopicAssignment;
	
	public int partNums = 0;
													
	public FeatureSequence wordSequence = null; // sequence representation of words in the document
	public FeatureSequence citationSequence = null;
	//add variables for different parts
	public FeatureSequence[] citationSequenceByPart = null;
												
	public TIntIntHashMap topic_counts = null; // bag representation of topic assignments
												
	public TIntIntHashMap word_counts = null; // bag representation of words in document

	public TIntIntHashMap links = null;
	//add variables for different parts
	public TIntIntHashMap[] linksByPart = null;
	
	public Vector<Integer> linksTopics = null;
	public String docName;
	public int docId;
	public int docLabel;
	public int docLength;

	public Document() {
		partNums = 2;
	}

	public Document(LabelSequence ts, FeatureSequence fs) {
		// assuming the content has been processed
		topicAssignments = ts;
		wordSequence = fs;
	}
	
	public void ProcessCSV(String name, int label, int docid, Vector<String> content,
			Vector<Citation> citations) {

		// create document object from content
		docName = name;
		docId = docid;
		docLabel = label;
		docLength = content.size();
		wordSequence = new FeatureSequence(Corpus.vocabulary);
		citationSequence = new FeatureSequence(Corpus.citationAlphabet);
		citationSequenceByPart = new FeatureSequence[partNums];
		for(int i = 0; i < citationSequenceByPart.length; i++){
			citationSequenceByPart[i] = new FeatureSequence(Corpus.citationAlphabet);
		}
		word_counts = new TIntIntHashMap();
		links = new TIntIntHashMap();
		linksByPart = new TIntIntHashMap[partNums];
		for(int i = 0; i < linksByPart.length; i++){
			linksByPart[i] = new TIntIntHashMap();
		}

		if (content == null) {
			//System.out.println("INFO: No content");
		} else {
			for (int i = 0; i < content.size(); i++) {
				String word = content.get(i);
				int index = Corpus.vocabulary.lookupIndex(word, true);
				wordSequence.add(index);
			}
			int[] words = wordSequence.getFeatures();
			for (int i = 0; i < wordSequence.length; i++) {
				word_counts.adjustOrPutValue(words[i], 1, 1);
			}
		}
		if (citations == null) {
			// System.out.println("INFO: No citations");
		} else {
			for (int i = 0; i < citations.size(); i++) {
				String citation = citations.get(i).id;
				int part = citations.get(i).part;
				int index = Corpus.citationAlphabet.lookupIndex(citation,true);
				citationSequence.add(index);
				citationSequenceByPart[part].add(index);
				links.adjustOrPutValue(index, 1, 1);
				linksByPart[part].adjustOrPutValue(index, 1, 1);
			}
		}
	}

	public void ProcessCSV_citaiton(String name, int label, int docid,
			Vector<String> citations) {

		// create document object from content
		docName = name;
		docId = docid;
		docLabel = label;
		
		links = new TIntIntHashMap();

		// String[] words = content.split(" ");
		// String[] edges = citations.split(" ");
		if (citations == null) {
			// System.out.println("INFO: No citations");

		} else {
			for (int i = 0; i < citations.size(); i++) {
				String edge = citations.get(i);
				int index = Corpus.citationAlphabet.lookupIndex(edge,true);
				//links.add(index);
				links.adjustOrPutValue(index, 1, 1);

			}
		}
	}

	private static final long serialVersionUID = 1;
	private static final int CURRENT_SERIAL_VERSION = 0;

}


