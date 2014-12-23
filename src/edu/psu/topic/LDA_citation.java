package edu.psu.topic;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Formatter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.TreeSet;
import java.util.Vector;
import java.util.regex.Pattern;
import java.io.*;

import gnu.trove.*;
import edu.psu.types.*;
import edu.psu.util.*;


public class LDA_citation extends Model{
	/**
	 * @param args
	 */
	public int numSamples = 0;
	//private static int numTopics;
	static ArrayList<Paper> paperList=new ArrayList<Paper>();
	static ArrayList<Paper> fullpaperList=new ArrayList<Paper>();

	private static LabelAlphabet newLabelAlphabet(int numTopics) {
		LabelAlphabet ret = new LabelAlphabet();
		for (int i = 0; i < numTopics; i++)
			ret.lookupIndex(i);
		return ret;
	}

	public Corpus readData(String Docdirectory,String citationDir) throws IOException, FileNotFoundException {
		InputStreamReader read;
		read = new InputStreamReader(new FileInputStream(new File(citationDir)),"UTF-8");	
        BufferedReader bufferedReader = new BufferedReader(read);
		int lineNum=0;
		String lineTxt = null;
		ArrayList<String> Docs=new ArrayList<String>();
		ArrayList<String> citations=new ArrayList<String>();
		
		while((lineTxt = bufferedReader.readLine()) != null)
		{
			lineNum++;
			if(lineNum%2==1)
			{Docs.add(lineTxt);
			}
			else
			{citations.add(lineTxt);}
		}
		
		bufferedReader.close();
		
		BufferedReader bufferedReader1 = new BufferedReader(new InputStreamReader(new FileInputStream(new File(Docdirectory)),"UTF-8"));
		String linePaper=null;
		 
		 
		while((linePaper=bufferedReader1.readLine())!=null)
		{
			String [] lines=linePaper.split("\\t");
			
			
			Paper paper=new Paper();
			paper.setID(lines[0]);
			paper.setTitle(lines[1]);
			paper.setYear(Integer.parseInt(lines[2]));
			fullpaperList.add(paper);
			if(Docs.contains(lines[0]))
			{
				//paper.setAbs(lines[3]);
				paperList.add(paper);
			}
		}
		
		bufferedReader1.close();
		
		Corpus corpus = new Corpus(paperList.size());
			
		for (int i = 0; i < paperList.size(); i++) {

			if (!Corpus.docAlphabet.contains(paperList.get(i).getID())) {
				Corpus.docAlphabet.lookupIndex(paperList.get(i).getID());
			}
		}			
		
		for (int i = 0; i < paperList.size(); i++) {

			ContextDocument doc = new ContextDocument();
			//doc.add_milestone(paperList.get(i).getID(), Corpus.docAlphabet.lookupIndex(paperList.get(i).getID()),(paperList.get(i).getTitle()+" "+paperList.get(i).getAbs()).toLowerCase(),Docs,citations);
			doc.add_milestone(paperList.get(i).getID(), Corpus.docAlphabet.lookupIndex(paperList.get(i).getID()),(paperList.get(i).getTitle()).toLowerCase(),Docs,citations);
			corpus.addDocument(doc,Corpus.docAlphabet.lookupIndex(paperList.get(i).getID()));
		}
		
			System.out.println("Total Cited Documents:"+ Corpus.citationAlphabet.size());
			System.out.println("Total Documents:" + corpus.docs.size());
			System.out.println("Total Vcabulary Size:"+ Corpus.vocabulary.size());
		return corpus;

	}

	public void InitializeParameters(Corpus corpus) {
		//corpus.numTypes = Corpus.vocabulary.size();

		//corpus.beta = 0.01;
		//corpus.betaSum = corpus.beta * corpus.numTypes;
		corpus.numCitations = Corpus.citationAlphabet.size();
		corpus.alpha = new double[numTopics];
		corpus.alphaSum = 0;
		for (int i = 0; i < numTopics; i++) {
			corpus.alpha[i] = 50.0 / (double) numTopics;
			corpus.alphaSum += corpus.alpha[i];
		}
		corpus.gamma = 0.01;
		corpus.gammaSum = corpus.gamma * corpus.numCitations;
		//corpus.typeTopicCounts = new TIntIntHashMap[corpus.numTypes];
		//corpus.tokensPerTopic = new int[numTopics];
		corpus.citationsPerTopic = new int[numTopics];
		//corpus.phi_train = new double[Corpus.vocabulary.size()][numTopics];
		corpus.theta_train = new double[Corpus.docAlphabet.size()][numTopics];
		corpus.psi_train = new double[Corpus.citationAlphabet.size()][numTopics];
		corpus.citationTopicCounts = new TIntIntHashMap[corpus.numCitations];
		for (int i = 0; i < corpus.citationTopicCounts.length; i++)
			corpus.citationTopicCounts[i] = new TIntIntHashMap();
		/*for (int i = 0; i < corpus.typeTopicCounts.length; i++)
			corpus.typeTopicCounts[i] = new TIntIntHashMap();*/

	}

	public void InitializeAssignments(Corpus corpus, LabelAlphabet topicAlphabet) {
		// initilize word, and citation factors
		this.random = new Randoms();

		// TODO AKM July 18: Why wasn't the next line there previously?
		// Only initialize the citation_topic_assignment
		
		for (int doc_index = 0; doc_index < corpus.docs.size(); doc_index++) {
			ContextDocument doc = corpus.docs.get(doc_index);
			//doc.citationTopicAssignment = new TIntIntHashMap();
			doc.citationTopicAssignmentDuplicate = new int[doc.citationSequence.getLength()];
			int[] keys = doc.citationSequence.getFeatures();
			for (int t = 0; t < doc.citationSequence.getLength(); t++) {
				int topic = r.nextInt(numTopics);
				//doc.citationTopicAssignment.put(keys[t], topic);// put docid and topic index  in citationtopicassignment
				doc.citationTopicAssignmentDuplicate[t]=topic;
				corpus.citationTopicCounts[keys[t]].adjustOrPutValue(topic, 1,1);
				corpus.citationsPerTopic[topic]++;
			}
		}

	}

	public void sampleOneDocument(Corpus corpus, ContextDocument doc) {
		TIntIntHashMap doc_topics = new TIntIntHashMap();

		int[] keys = doc.citationSequence.getFeatures();
		
		for (int t = 0; t < doc.citationSequence.getLength(); t++) 
		{
			//doc_topics.adjustOrPutValue(doc.citationTopicAssignment.get(keys[t]), 1, 1);
			doc_topics.adjustOrPutValue(doc.citationTopicAssignmentDuplicate[t], 1, 1);
		}
		keys = doc.citationSequence.getFeatures();
		
		for (int c = 0; c < doc.citationSequence.getLength(); c++) {
			//int topic = doc.citationTopicAssignment.get(keys[c]);
			int topic = doc.citationTopicAssignmentDuplicate[c];
			doc_topics.adjustOrPutValue(topic, -1, 0);
			corpus.citationTopicCounts[keys[c]].adjustOrPutValue(topic, -1, 0);
			corpus.citationsPerTopic[topic]--;
			TIntIntHashMap currentTypeTopicCounts = corpus.citationTopicCounts[keys[c]];
			double[] topicDistribution = new double[numTopics];
			double topicDistributionSum = 0, weight = 0;
			for (int t = 0; t < numTopics; t++) {
				weight = ((currentTypeTopicCounts.get(t) + corpus.gamma) / (corpus.citationsPerTopic[t] + corpus.gammaSum))
						* ((doc_topics.get(t) + corpus.alpha[t])/((doc.citationSequence.getLength()-1+corpus.alphaSum)));
				topicDistributionSum += weight;
				topicDistribution[t] = weight;
			}
			
			topic = random.nextDiscrete(topicDistribution, topicDistributionSum);
			//doc.citationTopicAssignment.put(keys[c], topic);
			doc.citationTopicAssignmentDuplicate[c] = topic;
			corpus.citationTopicCounts[keys[c]].adjustOrPutValue(topic, 1, 1);
			doc_topics.adjustOrPutValue(topic, 1, 1);
			corpus.citationsPerTopic[topic]++;
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

	public double testPerplexity(Corpus corpus) {
		double ll = 0;
		int sampleSize = 0;
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			ContextDocument document = (ContextDocument) corpus.docs.get(doc);
			sampleSize += document.citationSequence.getLength();
			
			double tmp_ll = 0;
			
			// For citation part
			int[] citations = document.citationSequence.getFeatures();
			for (int j = 0; j < citations.length; j++){
				for(int t = 0; t < numTopics; t++){
					tmp_ll += corpus.theta_train[doc][t]
							* corpus.psi_train[citations[j]][t];
					//System.out.println("theta_train : "+corpus.theta_train[doc][t]+", psi_train : "+corpus.psi_train[citations[j]][t]);
				}
				ll += Math.log(tmp_ll);
			}
		}
		//System.out.println("ll : "+ll+", sampleSize : "+sampleSize);
		return Math.exp(-1 * ll / sampleSize);
	}
	
	public void estimateParameters_single(Corpus corpus) {
		for(int doc = 0; doc < corpus.docs.size();doc++){
			// Initializing theta
			for(int i = 0;i < numTopics; i++) {
				corpus.theta_train[doc][i] = corpus.alpha[i]/(corpus.docs.get(doc).citationSequence.getLength()
						+corpus.alphaSum);
			}
			// Updating theta
			int[] topics = corpus.docs.get(doc).citationTopicAssignmentDuplicate;
			for(int i = 0;i < topics.length; i++){
				corpus.theta_train[doc][topics[i]] += 1.0/(corpus.docs.get(doc).citationSequence.getLength()
						+corpus.alphaSum);
			}
		}
		
		for (int k = 0; k < corpus.citationAlphabet.size(); k++) {
			for (int i = 0; i < numTopics; i++) {
				
				corpus.psi_train[k][i] = ((corpus.citationTopicCounts[k].get(i) + corpus.gamma) / (corpus.citationsPerTopic[i] + corpus.gammaSum));

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

	/*public static void main(String[] args) throws IOException{
		//String input = "D:/WorkPaceEclipse_2/citation_LDA-master/dataset/VLDB_TS_Occur/paper_info.txt";
		//String citationInput="D:/WorkPaceEclipse_2/citation_LDA-master/dataset/VLDB_TS_Occur/paper_citation.txt";
		String input = "D:/WorkPaceEclipse_2/citation_LDA-master/dataset/ACL/paper_ids-2010.txt";
		String citationInput="D:/WorkPaceEclipse_2/citation_LDA-master/dataset/ACL/Citing-Cited-2010.txt";
		int iterations = 5000;
		LDA_citation citationLda = new LDA_citation();
		System.out.println("Reading Data.....");
		Corpus corpus = citationLda.readData(input,citationInput);
		System.out.println("Done");
		//numTopics = 50;
		for(int i=0 ;i<=5;i++){	
			citationLda.numTopics = 50;
			citationLda.InitializeParameters(corpus);
			LabelAlphabet labelAlp=citationLda.newLabelAlphabet(citationLda.numTopics);
			
			citationLda.InitializeAssignments(corpus, labelAlp);
			citationLda.sampleCorpus(corpus, iterations, true);
			citationLda.estimateParameters_single(corpus);
		
			System.out.println(iterations+" : "+citationLda.testPerplexity(corpus));
		}
	}
	*/
	public static void main(String[] args) throws FileNotFoundException, IOException {
		// TODO Auto-generated method stub
		// Read documents
		String input = "D:/WorkPaceEclipse_2/citation_LDA-master/dataset/VLDB_TS/paper_info.txt";
		String citationInput="D:/WorkPaceEclipse_2/citation_LDA-master/dataset/VLDB_TS/paper_citation.txt";
		String output = "D:/WorkPaceEclipse_2/citation_LDA-master/dataset/VLDB_TS/result_citation_gamma_0.01_iter10000_#9.txt";
		
		int iterations = 10000;
		LDA_citation linkLda = new LDA_citation();
		System.out.println("Reading Data.....");
		Corpus corpus = linkLda.readData(input,citationInput);
		System.out.println("Done");
		//numTopics = 50;
		linkLda.numTopics = 50;
		linkLda.InitializeParameters(corpus);
	
		LabelAlphabet labelAlp=LDA_citation.newLabelAlphabet(linkLda.numTopics);
		linkLda.InitializeAssignments(corpus,labelAlp);
		linkLda.sampleCorpus(corpus, iterations, true);
		linkLda.estimateParameters_single(corpus);
		
		//ArrayList<TreeSet<IDSorter>> topicSortedWords = linkLda.getSortedWords(corpus,linkLda.numTopics);
		System.out.println("================================");
		Formatter out = new Formatter(new StringBuilder(), Locale.US);
		BufferedWriter bw1= new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(output),true)));
		
		bw1.write("# Topic_Citation");  
		bw1.newLine();
		HashMap<Integer,String> citationMap=new HashMap<Integer,String>();
		ArrayList<TreeSet<IDSorter>> topicSortedCitation = linkLda.getSortedCitation(corpus,linkLda.numTopics,citationMap);
		//double sum=0;
	
		for (int topic = 0; topic < linkLda.numTopics; topic++) {
			Iterator<IDSorter> iterator = topicSortedCitation.get(topic).iterator();
			
			out = new Formatter(new StringBuilder(), Locale.US);
			out.format ("%d\n",topic);
			int rank = 0;
			while (iterator.hasNext()&&rank<=20) {
				IDSorter idCountPair = iterator.next();
				//if(idCountPair.getWeight()>0.0001)
				//{
				String title=null;
				int year=-1;
				for(int i=0;i<fullpaperList.size();i++)
				{
					if(fullpaperList.get(i).getID().equals(citationMap.get(idCountPair.getID())))
					{
						title=fullpaperList.get(i).getTitle();
						year=fullpaperList.get(i).getYear();
					break;
					}
				}
				out.format("%s\t(%.4f)\t%s\t%d\n", citationMap.get(idCountPair.getID()), idCountPair.getWeight(),title,year);
				rank++;
			//}
			}
		//	System.out.println(out);
			String line=out.toString();
			bw1.write(line);  
			bw1.newLine();
		}
		bw1.write("# Doc_Topic");  
		bw1.newLine();
		for(int i = 0; i < corpus.docAlphabet.size(); i++){
			out = new Formatter(new StringBuilder(), Locale.US);
			String id = (String) corpus.docAlphabet.lookupObject(i);
			String title=null;
			int year=-1;
			for(int z=0;z<fullpaperList.size();z++)
			{
				if(fullpaperList.get(z).getID().equals(id))
				{
					title=fullpaperList.get(z).getTitle();
					year=fullpaperList.get(z).getYear();
					break;
				}
			}
			out.format("%s\t", corpus.docAlphabet.lookupObject(i));
			double max =0.0; int assign=-1;
			for(int topic = 0; topic<linkLda.numTopics; topic++){
				if(corpus.theta_train[i][topic]>max){
					max = corpus.theta_train[i][topic];
					assign=topic;
				}
				out.format ("%.4f\t",corpus.theta_train[i][topic]);
			}
			out.format("%d\t%.4f\t%s\t%d", assign,max,title,year);
			String line = out.toString();
			System.out.println(line);
			bw1.write(line);
			bw1.newLine();
		}
		out.close();
		bw1.close();

		System.out.println(iterations+" : "+linkLda.testPerplexity(corpus));
	}

}
