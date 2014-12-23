package edu.psu.topic;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Formatter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Locale;
import java.util.TreeSet;

import edu.psu.types.Alphabet;
import edu.psu.types.ContextDocument;
import edu.psu.types.Corpus;
import edu.psu.types.Dirichlet;
import edu.psu.types.FeatureSequence;
import edu.psu.types.IDSorter;
import edu.psu.types.LabelAlphabet;
import edu.psu.types.Paper;
import edu.psu.util.Randoms;
import gnu.trove.TIntIntHashMap;

public class collective_topic_model extends Model{
	/**
	 * @param args
	 */
	public int numSamples = 0;
	private static int numTopics;
	static int numCitations;
	public static Alphabet citationAlphabet = new Alphabet();
	public double gamma;
	public double gammaSum;
	public int[] citationsPerTopic;
	public double[][] psi_train;
	public TIntIntHashMap[] citationTopicCounts;
	static ArrayList<Paper> paperList=new ArrayList<Paper>();
	static ArrayList<Paper> fullpaperList=new ArrayList<Paper>();

	private static LabelAlphabet newLabelAlphabet(int numTopics) {
		LabelAlphabet ret = new LabelAlphabet();
		for (int i = 0; i < numTopics; i++)
			ret.lookupIndex(i);
		return ret;
	}

	public void readPaperInfor(String Docdirectory) throws IOException, FileNotFoundException 
	{
		 BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(Docdirectory)),"UTF-8"));
		 String linePaper=null;
		 
		 
		while((linePaper=bufferedReader.readLine())!=null)
		{
			String [] lines=linePaper.split("\\t");	
			Paper paper=new Paper();
			paper.setID(lines[0]);
			paper.setTitle(lines[1]);
			paper.setYear(Integer.parseInt(lines[2]));
			//System.out.println(paper.toString());
			fullpaperList.add(paper);
		}	
	}
	
	public Corpus readData(String citationDir) throws IOException, FileNotFoundException {
		 InputStreamReader read;
			read = new InputStreamReader(new FileInputStream(new File(citationDir)),"UTF-8");	
             BufferedReader bufferedReader = new BufferedReader(read);
		int lineNum=0;
		String lineTxt = null;
		ArrayList<String> Docs=new ArrayList();
		ArrayList<String> citations=new ArrayList();
		
		while((lineTxt = bufferedReader.readLine()) != null)
		{
			lineNum++;
			if(lineNum%2==1)
			{Docs.add(lineTxt);
			}
			else
			{citations.add(lineTxt);}
		}
		
		Corpus corpus = new Corpus(Docs.size());
		System.out.println(Docs.size()+" "+corpus.size);
		
		for (int i = 0; i < Docs.size(); i++) {

			if (!corpus.docNameAlphabet.contains(Docs.get(i))) {
				corpus.docNameAlphabet.lookupIndex(Docs.get(i));
			}
		}			
		
		for (int i = 0; i <  Docs.size(); i++) 
		{

			ContextDocument doc = new ContextDocument();
			doc.add_Auhtor(Docs.get(i), corpus.docNameAlphabet.lookupIndex(Docs.get(i)),citations.get(i));
			corpus.addDocument(doc,corpus.docNameAlphabet.lookupIndex(Docs.get(i)));
		}	
		
			System.out.println("Total  Documents:" + corpus.docs.size());
			
			return corpus;
	}

	public void InitializeParameters(Corpus corpus) {
		
		corpus.numTypes = corpus.vocabulary.size();
		corpus.beta = 0.01;
		corpus.betaSum = corpus.beta * corpus.numTypes;
		corpus.alpha = new double[numTopics];
		corpus.alphaSum = 0;
		for (int i = 0; i < numTopics; i++) {
			corpus.alpha[i] = 50.0 / (double) numTopics;
			corpus.alphaSum += corpus.alpha[i];
		}
		
		corpus.typeTopicCounts = new TIntIntHashMap[corpus.numTypes];
		corpus.tokensPerTopic = new int[numTopics];
		
		corpus.phi_train = new double[corpus.vocabulary.size()][numTopics];
		corpus.theta_train = new double[corpus.docAlphabet.size()][numTopics];
		
		for (int i = 0; i < corpus.typeTopicCounts.length; i++)
			corpus.typeTopicCounts[i] = new TIntIntHashMap();

	}

	//Initialize Author-Document and Venue-Document
	public void InitializeParameter_Author(Corpus corpus) {
		
		corpus.alpha = new double[numTopics];
		corpus.alphaSum = 0;
		for (int i = 0; i < numTopics; i++) {
			corpus.alpha[i] = 50.0 / (double) numTopics;
			corpus.alphaSum += corpus.alpha[i];
		}
		corpus.theta_train = new double[corpus.docAlphabet.size()][numTopics];
	}
	
	public void InitializeParameterAll(Corpus DocCorpus,Corpus AuthCorpus,Corpus VenueCorpus) 
	{
		InitializeParameter_Author(DocCorpus);
		InitializeParameter_Author(AuthCorpus);
		InitializeParameter_Author(VenueCorpus);
		numCitations = Corpus.citationAlphabet.size();
		gamma = 0.01;
		gammaSum = gamma * numCitations;
		citationsPerTopic = new int[numTopics];
		psi_train = new double[numCitations][numTopics];
		citationTopicCounts = new TIntIntHashMap[numCitations];
		for (int i = 0; i < citationTopicCounts.length; i++)
			citationTopicCounts[i] = new TIntIntHashMap();
		
	}
	
	public void InitializeAssignments(Corpus corpus, LabelAlphabet topicAlphabet) {
		// initilize word, and citation factors
		this.random = new Randoms();

		// TODO AKM July 18: Why wasn't the next line there previously?
		// this.typeTopicCounts = newTypeTopicCounts;
	
		for(int index=0;index<corpus.docs.size();index++)
		{
			ContextDocument doc = corpus.docs.get(index);
		 doc.citationTopicAssignment = new TIntIntHashMap();
		 int[] keys = doc.citationSet.keys();
		 for (int t = 0; t < keys.length; t++) {
			 int topic = r.nextInt(numTopics);
			doc.citationTopicAssignment.put(keys[t], topic);// put docid and topic index  in citationtopicassignment
			citationTopicCounts[keys[t]].adjustOrPutValue(topic, 1,1);
		     citationsPerTopic[topic]++;
		}
		}
	}

	public void sampleOneDocument(Corpus corpus, ContextDocument doc) {

		// decrement current sampling word

		// calculate document factor
		TIntIntHashMap doc_topics = new TIntIntHashMap();
		// sample for each position
		
		int[] topics = doc.wordTopicAssignment;
		int[] words = doc.wordSequence.getFeatures();
		//System.out.println(topics.length+" "+words.length+" "+doc.wordSequence.getLength());
		for (int i = 0; i < topics.length; i++) {
			doc_topics.adjustOrPutValue(doc.wordTopicAssignment[i], 1, 1);
		}
		//doc.citationTopicAssignment = new TIntIntHashMap();
		int[] keys = doc.citationTopicAssignment.keys();
		for (int t = 0; t < keys.length; t++) 
		{
			doc_topics.adjustOrPutValue(doc.citationTopicAssignment.get(keys[t]), 1, 1);
		}

		for (int j = 0; j <doc.wordSequence.getLength(); j++) {

				doc_topics.adjustOrPutValue(doc.wordTopicAssignment[j], -1, 0);
				corpus.typeTopicCounts[words[j]].adjustOrPutValue(doc.wordTopicAssignment[j],-1, 0);
				corpus.tokensPerTopic[doc.wordTopicAssignment[j]]--;
				// doc_topics.adjustOrPutValue(topics[j], -1, 0);
				TIntIntHashMap currentTypeTopicCounts = corpus.typeTopicCounts[words[j]];
				double[] topicDistribution = new double[numTopics];
				double topicDistributionSum = 0, weight = 0;
				for (int t = 0; t < numTopics; t++) {
					weight = ((currentTypeTopicCounts.get(t) + corpus.beta) / (corpus.tokensPerTopic[t] + corpus.betaSum))
							* ((doc_topics.get(t) + corpus.alpha[t]));
					topicDistributionSum += weight;
					topicDistribution[t] = weight;
				}
				// System.out.print("Sampled:");
				// System.out.println(random.nextDiscrete(topicDistribution,
				// topicDistributionSum));
				topics[j] = random.nextDiscrete(topicDistribution,topicDistributionSum);
				
				// topics[j] = r.nextInt(numTopics);
				doc.wordTopicAssignment[j] =topics[j];
				corpus.typeTopicCounts[words[j]].adjustOrPutValue(topics[j], 1,1);
				doc_topics.adjustOrPutValue(topics[j], 1, 1);
				corpus.tokensPerTopic[topics[j]]++;
		}
		keys = doc.citationTopicAssignment.keys();
		for (int c = 0; c < keys.length; c++) {
			int topic = doc.citationTopicAssignment.get(keys[c]);
			doc_topics.adjustOrPutValue(topic, -1, 0);
			citationTopicCounts[keys[c]].adjustOrPutValue(topic, -1, 0);
			citationsPerTopic[topic]--;
			TIntIntHashMap currentTypeTopicCounts = citationTopicCounts[keys[c]];
			double[] topicDistribution = new double[numTopics];
			double topicDistributionSum = 0, weight = 0;
			for (int t = 0; t < numTopics; t++) {
				weight = ((currentTypeTopicCounts.get(t) + gamma) / (citationsPerTopic[t] + gammaSum))
						* ((doc_topics.get(t) + corpus.alpha[t]));
				topicDistributionSum += weight;
				topicDistribution[t] = weight;
			}
			// System.out.print("Sampled:");
			// System.out.println(random.nextDiscrete(topicDistribution,
			// topicDistributionSum));
			topic = random.nextDiscrete(topicDistribution, topicDistributionSum);
			// topics[j] = r.nextInt(numTopics);
			doc.citationTopicAssignment.put(keys[c], topic);
			citationTopicCounts[keys[c]].adjustOrPutValue(topic, 1, 1);
			doc_topics.adjustOrPutValue(topic, 1, 1);
			citationsPerTopic[topic]++;
		}

	}

	public void sampleOneDocument_author(Corpus corpus, ContextDocument doc) {

		// decrement current sampling word

		// calculate document factor
		TIntIntHashMap doc_topics = new TIntIntHashMap();
		// sample for each position
		
		//doc.citationTopicAssignment = new TIntIntHashMap();
		int[] keys = doc.citationTopicAssignment.keys();
		System.out.println(doc.docName);
		for (int t = 0; t < keys.length; t++) 
		{
			doc_topics.adjustOrPutValue(doc.citationTopicAssignment.get(keys[t]), 1, 1);
		}

		keys = doc.citationTopicAssignment.keys();
		for (int c = 0; c < keys.length; c++) {
			int topic = doc.citationTopicAssignment.get(keys[c]);
			doc_topics.adjustOrPutValue(topic, -1, 0);
			citationTopicCounts[keys[c]].adjustOrPutValue(topic, -1, 0);
			citationsPerTopic[topic]--;
			TIntIntHashMap currentTypeTopicCounts = citationTopicCounts[keys[c]];
			double[] topicDistribution = new double[numTopics];
			double topicDistributionSum = 0, weight = 0;
			for (int t = 0; t < numTopics; t++) {
				weight = ((currentTypeTopicCounts.get(t) + gamma) / (citationsPerTopic[t] + gammaSum))
						* ((doc_topics.get(t) + corpus.alpha[t]));
				topicDistributionSum += weight;
				topicDistribution[t] = weight;
			}
			// System.out.print("Sampled:");
			// System.out.println(random.nextDiscrete(topicDistribution,
			// topicDistributionSum));
			topic = random.nextDiscrete(topicDistribution, topicDistributionSum);
			// topics[j] = r.nextInt(numTopics);
			doc.citationTopicAssignment.put(keys[c], topic);
			citationTopicCounts[keys[c]].adjustOrPutValue(topic, 1, 1);
			doc_topics.adjustOrPutValue(topic, 1, 1);
			citationsPerTopic[topic]++;
		}

	}
	public void sampleColletiveCorpus(Corpus DocCorpus,Corpus AuthCorpus,Corpus VenueCorpus,int iterations)
	{
		for(int ite=0;ite<iterations;ite++)
		{
		
		for(int j=0;j<AuthCorpus.docs.size();j++)
		{
			sampleOneDocument_author(AuthCorpus,(ContextDocument) AuthCorpus.docs.get(j));
		}

		for(int k=0;k<VenueCorpus.docs.size();k++)
		{
			sampleOneDocument_author(VenueCorpus,(ContextDocument) VenueCorpus.docs.get(k));
		}
		for (int i = 0; i < DocCorpus.docs.size(); i++) {
			sampleOneDocument_author(DocCorpus,(ContextDocument) DocCorpus.docs.get(i));
		}
		
		}
		
	}
	
	public void sampleCorpus(Corpus corpus, int iterations,boolean isContextAware) {
		for(int ite=0;ite<iterations;ite++)
		{
		for (int i = 0; i < corpus.docs.size(); i++) {
			sampleOneDocument(corpus,(ContextDocument) corpus.docs.get(i));
		}
		}
	}


	public double sampleLikelihood(int numSamples, Corpus corpus) {
		double ll = 0;
		for (int doc = 0; doc < corpus.docs.size(); doc++) {

			TIntIntHashMap word_counts = corpus.docs.get(doc).word_counts;
			int[] keys = word_counts.keys();
			double tmp_ll = 0;
			for (int i = 0; i < keys.length; i++) {
				for (int t = 0; t < numTopics; t++) {
					tmp_ll += corpus.theta_train[doc][t]
							* corpus.phi_train[keys[i]][t];
				}
				ll += Math.log(tmp_ll * word_counts.get(keys[i]));
			}

		}
		return ll;
	}

	public double sampleLikelihood(int numSamples, Corpus corpus,
			TIntIntHashMap split) {
		double ll = 0;
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			if (!split.contains(doc)) {
				TIntIntHashMap word_counts = corpus.docs.get(doc).word_counts;
				int[] keys = word_counts.keys();
				double tmp_ll = 0;
				for (int i = 0; i < keys.length; i++) {
					for (int t = 0; t < numTopics; t++) {
						tmp_ll += corpus.theta_train[doc][t]
								* corpus.phi_train[keys[i]][t];
					}
					ll += Math.log(tmp_ll * word_counts.get(keys[i]));
				}
			}
		}
		return ll;
	}

	public double empiricalLikelihood(int numSamples, Corpus corpus) {
		double[][] likelihoods = new double[corpus.docs.size()][numSamples];
		double[] multinomial = new double[corpus.numTypes];
		double[] topicDistribution, currentSample, currentWeights;
		Dirichlet topicPrior = new Dirichlet(corpus.alpha);

		int sample, doc, topic, type, token, seqLen;
		FeatureSequence fs;

		for (sample = 0; sample < numSamples; sample++) {
			topicDistribution = topicPrior.nextDistribution();
			Arrays.fill(multinomial, 0.0);

			for (topic = 0; topic < numTopics; topic++) {
				for (type = 0; type < corpus.numTypes; type++) {
					multinomial[type] += topicDistribution[topic]
							* (corpus.beta + corpus.typeTopicCounts[type]
									.get(topic))
							/ (corpus.betaSum + corpus.tokensPerTopic[topic]);
				}
			}

			// Convert to log probabilities
			for (type = 0; type < corpus.numTypes; type++) {
				assert (multinomial[type] > 0.0);
				multinomial[type] = Math.log(multinomial[type]);
			}

			for (doc = 0; doc < corpus.docs.size(); doc++) {
				fs = (FeatureSequence) corpus.docs.get(doc).wordSequence;
				seqLen = fs.getLength();

				for (token = 0; token < seqLen; token++) {
					type = fs.getIndexAtPosition(token);

					// Adding this check since testing instances may
					// have types not found in training instances,
					// as pointed out by Steven Bethard.
					if (type < corpus.numTypes) {
						likelihoods[doc][sample] += multinomial[type];
					}
				}
			}
		}

		double averageLogLikelihood = 0.0;
		double logNumSamples = Math.log(numSamples);
		for (doc = 0; doc < corpus.docs.size(); doc++) {
			double max = Double.NEGATIVE_INFINITY;
			for (sample = 0; sample < numSamples; sample++) {
				if (likelihoods[doc][sample] > max) {
					max = likelihoods[doc][sample];
				}
			}

			double sum = 0.0;
			for (sample = 0; sample < numSamples; sample++) {
				sum += Math.exp(likelihoods[doc][sample] - max);
			}

			averageLogLikelihood += Math.log(sum) + max - logNumSamples;
		}

		return averageLogLikelihood;

	}

	public void precision(Corpus corpus) {

	}
	
public ArrayList<TreeSet<IDSorter>> getSortedCitation (int numTopics,HashMap<Integer,String> citationMap) {
		
		ArrayList<TreeSet<IDSorter>> topicSortedCitation = new ArrayList<TreeSet<IDSorter>>(numTopics);

		// Initialize the tree sets
		for (int topic = 0; topic < numTopics; topic++) {
			topicSortedCitation.add(new TreeSet<IDSorter>());
		}

		// Collect counts
		ArrayList<String> entries=Corpus.citationAlphabet.entries;
		for (int type = 0; type <Corpus.citationAlphabet.size(); type++) {
			int index=Corpus.citationAlphabet.lookupIndex(entries.get(type));
			citationMap.put(index, entries.get(type));
			for (int i = 0; i < numTopics; i++) {			
				psi_train[type][i] = ((citationTopicCounts[index].get(i) + gamma) / (citationsPerTopic[i] + gammaSum));
				topicSortedCitation.get(i).add(new IDSorter(index, psi_train[type][i]));
			}

		}


		return topicSortedCitation;
	}
	

	public static void main(String[] args) throws FileNotFoundException, IOException {
		// TODO Auto-generated method stub
		// Read documents
		String Docinput = "E:/project/citation_LDA-master/dataset/PLSA_dataset/paper_ids.txt";
		String Authinput="E:/project/citation_LDA-master/dataset/PLSA_dataset/author_ids.txt";
		String Venueinput="E:/project/citation_LDA-master/dataset/PLSA_dataset/VenueIDs.txt";
		String Doc_citation="E:/project/citation_LDA-master/dataset/PLSA_dataset/DocRefer.txt";
		String Auth_citation="E:/project/citation_LDA-master/dataset/PLSA_dataset/AuthorCiti.txt";
		String Venue_citation="E:/project/citation_LDA-master/dataset/PLSA_dataset/VenueRefer.txt";
		String output = "E:/project/Mallet/data/ACLresult/colletive_model_authFirst.txt";
		
		int iterations = 1000;
		collective_topic_model linkLda = new collective_topic_model();
		System.out.println("Reading Data.....");
		linkLda.readPaperInfor(Docinput);
		Corpus Doc=linkLda.readData(Doc_citation);
		Corpus Auth=linkLda.readData(Auth_citation);
		Corpus Venue=linkLda.readData(Venue_citation);
		
		System.out.println("Done");
		numTopics = 100;
		linkLda.InitializeParameterAll(Doc,Auth,Venue);
	
		LabelAlphabet labelAlp=linkLda.newLabelAlphabet(numTopics);
		linkLda.InitializeAssignments(Doc,labelAlp);
		linkLda.InitializeAssignments(Auth,labelAlp);
		linkLda.InitializeAssignments(Venue,labelAlp);
		
		linkLda.sampleColletiveCorpus(Doc,Auth,Venue,iterations);

		
	//	ArrayList<TreeSet<IDSorter>> topicSortedWords = linkLda.getSortedWords(DocCorpus,numTopics);
		System.out.println("================================");
		Formatter out = new Formatter(new StringBuilder(), Locale.US);
		BufferedWriter bw1= new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(output),true)));
		/*
		bw1.write("# Topic_word");  
		bw1.newLine();
		for (int topic = 0; topic < numTopics; topic++) {
			Iterator<IDSorter> iterator = topicSortedWords.get(topic).iterator();
			
			out = new Formatter(new StringBuilder(), Locale.US);
			out.format ("%d ",topic);
			int rank = 0;
			while (iterator.hasNext() && rank<20) {
				IDSorter idCountPair = iterator.next();
				out.format("%s (%.4f) ", DocCorpus.vocabulary.lookupObject(idCountPair.getID()), idCountPair.getWeight());
				rank++;
			}
		  System.out.println(out);
			String line=out.toString();
			bw1.write(line);  
			bw1.newLine();
		}
		*/
		bw1.write("# Topic_Citation");  
		bw1.newLine();
		HashMap<Integer,String> citationMap=new HashMap<Integer,String>();
		ArrayList<TreeSet<IDSorter>> topicSortedCitation = linkLda.getSortedCitation(numTopics,citationMap);
		double sum=0;
	
		for (int topic = 0; topic < numTopics; topic++) {
			Iterator<IDSorter> iterator = topicSortedCitation.get(topic).iterator();
			
			out = new Formatter(new StringBuilder(), Locale.US);
			out.format ("%d ",topic);
			int rank = 0;
			while (iterator.hasNext()&& rank<20) {
				IDSorter idCountPair = iterator.next();
				String title=null;
				for(int i=0;i<fullpaperList.size();i++)
				{
					if(fullpaperList.get(i).getID().equals(citationMap.get(idCountPair.getID())))
					{title=fullpaperList.get(i).getTitle();
					break;
					}
				}
				out.format("%s (%.4f) %s \n", citationMap.get(idCountPair.getID()), idCountPair.getWeight(),title);
				rank++;
			}
		//	System.out.println(out);
			String line=out.toString();
			bw1.write(line);  
			bw1.newLine();
		}	
		bw1.close();
	}
}
