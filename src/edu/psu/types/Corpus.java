package edu.psu.types;

import java.util.*;

import edu.psu.types.*;
import edu.psu.util.ReadDirectory;
import gnu.trove.*;

import java.io.*;
import java.math.*;

public class Corpus implements Serializable {
	public static Alphabet vocabulary = new Alphabet();
	public static Alphabet citationAlphabet = new Alphabet();
	public static Alphabet docAlphabet = new Alphabet();  //���浠ュ����炬��妗ｅ��瀛�
	public  Alphabet docNameAlphabet = new Alphabet();  
	public static ArrayList<String> citationVector=new ArrayList<String>();

	public static int size;
	public Vector<ContextDocument> docs;
	public TIntObjectHashMap<ContextDocument> docsTable;
	//public Vector<WordOccurrences> words;
	//public static TIntDoubleHashMap idf;
	public TIntIntHashMap[] typeTopicCounts;
	public TIntIntHashMap[] citationTopicCounts;
	
	public int numTypes;
	public int numCitations;
	//add numPart for representing the num of parts
	public int numParts;
	//public int[] tokensPerTopic;
	public int[] citationsPerTopic;
	//add citation topic for different parts
	public int[][] citationsPerTopicByPart;

	public Corpus() {

	}

	public Corpus(int capacity) {
		docs = new Vector<ContextDocument>();
		for (int i = 0; i < capacity; i++) {
			docs.add(new ContextDocument());
		}
		size = capacity;
	}

	public void addDocument(ContextDocument doc) {
		if (docs == null) {
			docs = new Vector<ContextDocument>();
		}
		docs.add(doc);
		size = docs.size();
	}

	public void addDocument(ContextDocument doc, int index) {
		docs.set(index, doc);
		if (docsTable == null) {
			docsTable = new TIntObjectHashMap<ContextDocument>();
		}
		docsTable.put(index, doc);
	}

	public Document getDocument(int index) {
		if (docs.size() <= index) {
			System.err.println("doc id exceed corpus size!");
			System.exit(-1);
		}
		return docs.get(index);
	}

	public void addDocuments(Vector<ContextDocument> documents) {

		docs = documents;
		size = docs.size();
	}

	private static final long serialVersionUID = 1;
	private static final int CURRENT_SERIAL_VERSION = 0;


}
