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
	//public static Alphabet citedAlphabet=new Alphabet();
	public static Alphabet docAlphabet = new Alphabet();  //可以存放文档名字
	public  Alphabet docNameAlphabet = new Alphabet();  
	public static ArrayList<String> citationVector=new ArrayList<String>();
	/*public TIntDoubleHashMap[] cluster_states;
	public TIntDoubleHashMap[] cluster_probs_reconst;
	public static int citationDocIndex=0;
	public static double[] mean;
	public double[][] eignVector;
	public double[][] normalizedEignVector;
	public double[][][] subspaceEignVector;
	public double[] mean_projection;
	public double[][] c_proj;
	public double[][] cluster;
	public int[] cluster_doc;
	public int[] cluster_size;
	public double[][] belong;
	public double[] cluster_probs;
	public TIntIntHashMap[] word_class_occ;*/

	//public double[] alphas;
	//public TIntDoubleHashMap[] subspace_vectors;
	//public int[] doc_cluster;
	//public int vocab_size;
	public static int size;
	public Vector<ContextDocument> docs;
	public TIntObjectHashMap<ContextDocument> docsTable;
	//public Vector<WordOccurrences> words;
	//public static TIntDoubleHashMap idf;
	public TIntIntHashMap[] typeTopicCounts;
	public TIntIntHashMap[] citationTopicCounts;
	
	public int numTypes;
	public int numCitations;
	public int[] tokensPerTopic;
	public int[] citationsPerTopic;
	
	//public boolean isPerplex = false;
	//public static TIntIntHashMap split;
	//public static Map<Integer, Vector<ContextDocument>> test_docs = new HashMap<Integer, Vector<ContextDocument>>();
	//public static Map<Integer, Vector<ContextDocument>> train_docs = new HashMap<Integer, Vector<ContextDocument>>();
	//public static Map<Integer, ConceptVectors> concepts = new HashMap<Integer, ConceptVectors>();
	//public TIntIntHashMap[] trainCitations;
	//public Vector<RankedList> citedTopicLists;
	public static RankedList citedList = new RankedList();
	public int maxTokens;

	public Corpus() {

	}

	public class Splits {
		TIntIntHashMap[] tests;
		TIntIntHashMap[] trains;

		public Splits(int numFolds) {
			tests = new TIntIntHashMap[numFolds];
			trains = new TIntIntHashMap[numFolds];
			for (int i = 0; i < numFolds; i++) {
				tests[i] = new TIntIntHashMap();
				trains[i] = new TIntIntHashMap();
			}
		}

		public void produceSplits(Integer[] instances, int numFolds) {
			// shuffle array
			Random r = new Random(21);
			System.out.println(instances.length);
			for (int i = 0; i < instances.length; i++) {
				int rand = r.nextInt(instances.length - i) + i;
				int temp = instances[i];
				instances[i] = instances[rand];
				instances[rand] = temp;
			}
			for (int i = 0; i < instances.length; i++) {
				int which_split = i % numFolds;
				// System.out.println(instances[i].intValue());

				tests[which_split].put(instances[i].intValue(), 1);
				for (int j = 0; j < numFolds; j++) {
					if (j != which_split) {
						trains[j].put(instances[i].intValue(), 1);
					}
				}

			}

		}

	}

	/*public void prepareSplits(int numFolds, boolean isCitationInformation) {

		Integer[] entries = new Integer[docs.size()];
		for (int i = 0; i < docs.size(); i++) {
			entries[i] = i;
		}
		trainCitations = new TIntIntHashMap[numFolds];
		for (int i = 0; i < numFolds; i++) {
			Vector<ContextDocument> test = new Vector<ContextDocument>();
			Vector<ContextDocument> train = new Vector<ContextDocument>();
			test_docs.put(i, test);
			train_docs.put(i, train);
		}
		Splits split = new Splits(numFolds);
		split.produceSplits(entries, numFolds);
		for (int i = 0; i < numFolds; i++) {
			trainCitations[i] = new TIntIntHashMap();
			int[] keys = split.trains[i].keys();
			for (int j = 0; j < keys.length; j++) {
				ContextDocument doc = (ContextDocument) this.docs.get(keys[j]);
				int[] cited = doc.citationSet.keys();
				for (int k = 0; k < cited.length; k++) {
					trainCitations[i].put(cited[k], 1);
				}

			}
			keys = split.tests[i].keys();
			System.out
					.println("Before removing doc with no cited articles in fold"
							+ i
							+ ": train="
							+ split.trains[i].size()
							+ ":Test=" + split.tests[i].size());
			// take documents away from test set which dont have citations.
			for (int j = 0; j < keys.length; j++) {
				ContextDocument doc = (ContextDocument) this.docs.get(keys[j]);
				if (doc.citationSet.size() == 0) {
					split.trains[i].put(keys[j], 1);
					split.tests[i].remove(keys[j]);
				} else {
					int[] cited = doc.citationSet.keys();
					boolean seen_in_training = false;
					for (int k = 0; k < cited.length; k++) {
						if (trainCitations[i].containsKey(cited[k])) {
							seen_in_training = true;
						}
					}
					if (!seen_in_training) {
						split.trains[i].put(keys[j], 1);
						split.tests[i].remove(keys[j]);
					}
				}
			}
			System.out
					.println("After removing doc with no cited articles in fold"
							+ i
							+ ": train="
							+ split.trains[i].size()
							+ ":Test=" + split.tests[i].size());

		}

		// remove all the citation from test which do not appear in training.

		for (int i = 0; i < numFolds; i++) {
			Vector<ContextDocument> test = test_docs.get(i);
			Vector<ContextDocument> train = train_docs.get(i);
			int[] Keys = split.tests[i].keys();
			for (int j = 0; j < Keys.length; j++) {
				ContextDocument doc = (ContextDocument) docs.get(Keys[j]);
				test.add(doc.RemoveCitation(trainCitations[i]));
			}
			Keys = split.trains[i].keys();
			for (int j = 0; j < Keys.length; j++) {
				train.add(docs.get(Keys[j]));
			}
			// remove citation information
			// clone the documents in test and train splits and
			// remove citation information from test documents
			// do for citation recommendation only--check to see if test
			// document has any outgoing citation link or not

		}
		if (!isCitationInformation) {
			for (int i = 0; i < numFolds; i++) {
				Vector<ContextDocument> test = test_docs.get(i);

				for (int j = 0; j < test.size(); j++) {
					// System.out.print(((ContextDocument)test.get(j)).docId+":"+((ContextDocument)test.get(j)).citationSet.size()+":");
					ContextDocument doc = test.get(j);

					doc = doc.RemoveCitation();

					test.remove(j);
					test.insertElementAt(doc, j);
					// System.out.println(((ContextDocument)test.get(j)).docId+":"+((ContextDocument)test.get(j)).citationSet.size()+":");
				}
				// remove citation information
				// clone the documents in test and train splits and
				// remove citation information from test documents
				// do for citation recommendation only--check to see if test
				// document has any outgoing citation link or not

			}
		}

	}*/

	/*public static Corpus assignDocumentIds(String directory) {
		String[] files = ReadDirectory.list(directory);
		System.out.println("There are total " + files.length + " Files in "
				+ directory);
		for (int i = 0; i < files.length; i++) {
			String name = files[i].split("/")[files[i].split("/").length - 1]
					.replace(".txt", "");// check!!!
			if (!Corpus.docAlphabet.contains(name)) {
				Corpus.docAlphabet.lookupIndex(name);
			}
		}
		Corpus corpus = new Corpus(Corpus.docAlphabet.size());
		return corpus;

	}
*/
	
	public void readData(String directory) {
		String[] files = ReadDirectory.list(directory);
		System.out.println("There are total " + files.length + " Files in "
				+ directory);
		String line = null;
		maxTokens=0;

		try {
			for (int i = 0; i < files.length; i++) {
				String name = files[i].split("/")[files[i].split("/").length - 1]
						.replace(".txt", "");// check!!!
				FileReader fr = new FileReader(files[i]);
				BufferedReader br = new BufferedReader(fr);
				try {
					Vector<String> lines = new Vector<String>();
					while ((line = br.readLine()) != null) {
						lines.add(line);
					}
					ContextDocument doc = new ContextDocument();
					doc.add(name, Corpus.docAlphabet.lookupIndex(name), lines);
					this.addDocument(doc, Corpus.docAlphabet.lookupIndex(name));
					if (maxTokens < doc.wordLength)
						maxTokens = doc.wordLength;
					fr.close();
					br.close();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}

			}
			for (int i = 0; i < docsTable.size(); i++) {
				ContextDocument doc = (ContextDocument) docsTable.get(i);
				for (int j = 0; j < doc.citationSet.size(); j++) {
					Corpus.citedList.add(doc.citationSet.get(j), 1);
				}
			}

			FileWriter fw = new FileWriter("./cited_docs");
			BufferedWriter bw = new BufferedWriter(fw);

			for (int i = 0; i < Corpus.citedList.indices.size(); i++) {
				bw.write(Corpus.docAlphabet
						.lookupObject(Corpus.citedList.indices.get(i))
						+ ":"
						+ Corpus.citedList.values.get(Corpus.citedList.indices
								.get(i)) + "\n");
			}
			bw.close();
			fw.close();

			System.out.println("Total Cited Documents:"
					+ Corpus.citationAlphabet.size());
			System.out.println("Total Documents:" + this.docs.size());
			System.out.println("Total Vcabulary Size:"
					+ Corpus.vocabulary.size());
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void setWindowLength(int length) {
		for (int i = 0; i < this.docs.size(); i++) {
			ContextDocument doc = (ContextDocument) this.docs.get(i);

			TIntObjectHashMap<Vector<Integer>> contextObject = doc.contextObject;
			int[] keys = contextObject.keys();
			for (int k = 0; k < keys.length; k++) {
				Vector<Integer> tmpContext = contextObject.get(keys[k]);
				for (int j = 1; j <= length; j++) {
					if (keys[k] + j < doc.sentenseLength) {
						contextObject.put(keys[k] + j,
								(Vector<Integer>) tmpContext.clone());
					}
					if (keys[k] - j >= 0) {
						contextObject.put(keys[k] - j,
								(Vector<Integer>) tmpContext.clone());
					}
				}
			}
		}
	}

	public void adaptiveWindowLength(int fold) {
		// assimilate topic assignment from present context, context+1,
		// context-1 and cited document's initial 10 sentences.
		System.out.print("Updating the context windows.......");
		int total = 0;
		double avgContextLength = 0, avgLeftContextLength = 0, avgRightContextLength = 0;
		for (int i = 0; i < this.docs.size(); i++) {
			ContextDocument doc = (ContextDocument) this.docs.get(i);
			doc.setWindowLenghtAdaptive(this, fold);
			if (doc.avgContextLength > 0) {
				total++;
			}
			avgContextLength += doc.avgContextLength;
			avgLeftContextLength += doc.avgLeftContextLength;
			avgRightContextLength += doc.avgRightContextLength;
		}
		avgContextLength /= total;
		avgLeftContextLength /= total;
		avgRightContextLength /= total;
		System.out.println("Done!:Avg Context Length=" + avgContextLength
				+ ":Avg Left Context Length=" + avgLeftContextLength
				+ ":Avg Right Context Length=" + avgRightContextLength);
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
