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


public class linkLDA_milestone_test extends Model{
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
			//System.out.println(paper.toString());
			fullpaperList.add(paper);
			if(Docs.contains(lines[0]))
			{
				paper.setAbs(lines[3]);
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
			doc.add_milestone(paperList.get(i).getID(), Corpus.docAlphabet.lookupIndex(paperList.get(i).getID()),(paperList.get(i).getTitle()+" "+paperList.get(i).getAbs()).toLowerCase(),Docs,citations);
			corpus.addDocument(doc,Corpus.docAlphabet.lookupIndex(paperList.get(i).getID()));
		}
		
			System.out.println("Total Cited Documents:"+ Corpus.citationAlphabet.size());
			System.out.println("Total Documents:" + corpus.docs.size());
			System.out.println("Total Vcabulary Size:"+ Corpus.vocabulary.size());
		return corpus;

	}

	public void InitializeParameters(Corpus corpus) {
		corpus.numTypes = corpus.vocabulary.size();

		corpus.beta = 0.01;
		corpus.betaSum = corpus.beta * corpus.numTypes;
		corpus.numCitations = corpus.citationAlphabet.size();
		corpus.alpha = new double[numTopics];
		corpus.alphaSum = 0;
		for (int i = 0; i < numTopics; i++) {
			corpus.alpha[i] = 50.0 / (double) numTopics;
			corpus.alphaSum += corpus.alpha[i];
		}
		corpus.gamma = 0.01;
		corpus.gammaSum = corpus.gamma * corpus.numCitations;
		corpus.typeTopicCounts = new TIntIntHashMap[corpus.numTypes];
		corpus.tokensPerTopic = new int[numTopics];
		corpus.citationsPerTopic = new int[numTopics];
		corpus.phi_train = new double[corpus.vocabulary.size()][numTopics];
		corpus.theta_train = new double[corpus.docAlphabet.size()][numTopics];
		corpus.psi_train = new double[corpus.citationAlphabet.size()][numTopics];
		//citationTopicCounts(大小不加docAlphabet应该也是对的)
		corpus.citationTopicCounts = new TIntIntHashMap[corpus.numCitations];
		for (int i = 0; i < corpus.citationTopicCounts.length; i++)
			corpus.citationTopicCounts[i] = new TIntIntHashMap();
		for (int i = 0; i < corpus.typeTopicCounts.length; i++)
			corpus.typeTopicCounts[i] = new TIntIntHashMap();

	}

	public void InitializeAssignments(Corpus corpus, LabelAlphabet topicAlphabet) {
		// initilize word, and citation factors
		this.random = new Randoms();

		// TODO AKM July 18: Why wasn't the next line there previously?
		// this.typeTopicCounts = newTypeTopicCounts;

		for (int doc_index = 0; doc_index < corpus.docs.size(); doc_index++) {
			//System.out.println("doc_index "+doc_index);
			ContextDocument doc = corpus.docs.get(doc_index);
			//System.out.println(doc.wordSequence.getLength());
			doc.wordTopicAssignment = new int[doc.wordSequence.getLength()];
			doc.citationTopicAssignmentDuplicate = new int[doc.citationSequence.getLength()];

			int[] words = doc.wordSequence.getFeatures();
			/*
			System.out.println(doc_index+" "+topics.length);
			for(int windex=0;windex<topics.length;windex++)
			{
				System.out.println(topics[windex]);
			}
			*/
			// put a check to see if wordlength exceeds words.length
			for (int i = 0; i < doc.wordSequence.getLength(); i++) {

					int topic = r.nextInt(numTopics);
					doc.wordTopicAssignment[i]=topic;
					
					corpus.typeTopicCounts[words[i]].adjustOrPutValue(topic, 1, 1);
					corpus.tokensPerTopic[topic]++;	
			}
			
			//only works with citation set(not consider occurrance)
			//for citation occurrance, NEED CHANGE THE CODE
			doc.citationTopicAssignment = new TIntIntHashMap();
			int[] keys = doc.citationSequence.getFeatures();
			//int[] keys = doc.citationSet.keys();
			for (int t = 0; t < doc.citationSequence.getLength(); t++) {
				int topic = r.nextInt(numTopics);
				doc.citationTopicAssignment.put(keys[t], topic);// put docid and topic index  in citationtopicassignment
				doc.citationTopicAssignmentDuplicate[t]=topic;
				corpus.citationTopicCounts[keys[t]].adjustOrPutValue(topic, 1,1);
				corpus.citationsPerTopic[topic]++;
			}
		}

	}

	public void sampleOneDocument(Corpus corpus, ContextDocument doc) {

		// decrement current sampling word

		// calculate document factor
		TIntIntHashMap doc_topics = new TIntIntHashMap();  //topic, and the occurrance time of the topic
		// sample for each position
		
		int[] topics = doc.wordTopicAssignment;
		int[] words = doc.wordSequence.getFeatures();
		//System.out.println(topics.length+" "+words.length+" "+doc.wordSequence.getLength());
		for (int i = 0; i < topics.length; i++) {
			doc_topics.adjustOrPutValue(doc.wordTopicAssignment[i], 1, 1);
		}
		//doc.citationTopicAssignment = new TIntIntHashMap();
		//int[] keys = doc.citationTopicAssignment.keys();
		int[] keys = doc.citationSequence.getFeatures();
		for (int t = 0; t < doc.citationSequence.getLength(); t++) 
		{
			//doc_topics.adjustOrPutValue(doc.citationTopicAssignment.get(keys[t]), 1, 1);
			doc_topics.adjustOrPutValue(doc.citationTopicAssignmentDuplicate[t], 2, 1);
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
							* ((doc_topics.get(t) + corpus.alpha[t])/(doc.docLength+doc.citationSequence.getLength()-1+corpus.alphaSum));
					topicDistributionSum += weight;
					topicDistribution[t] = weight;
				}
			
				topics[j] = random.nextDiscrete(topicDistribution,topicDistributionSum);
				
				doc.wordTopicAssignment[j] =topics[j];
				corpus.typeTopicCounts[words[j]].adjustOrPutValue(topics[j], 1,1);
				doc_topics.adjustOrPutValue(topics[j], 1, 1);
				corpus.tokensPerTopic[topics[j]]++;
		}
		//keys = doc.citationTopicAssignment.keys();
		keys = doc.citationSequence.getFeatures();
		for (int c = 0; c < doc.citationSequence.getLength(); c++) {
			//int topic = doc.citationTopicAssignment.get(keys[c]);
			int topic = doc.citationTopicAssignmentDuplicate[c];
			doc_topics.adjustOrPutValue(topic, -2, 0);
			corpus.citationTopicCounts[keys[c]].adjustOrPutValue(topic, -1, 0);
			corpus.citationsPerTopic[topic]--;
			TIntIntHashMap currentTypeTopicCounts = corpus.citationTopicCounts[keys[c]];
			double[] topicDistribution = new double[numTopics];
			double topicDistributionSum = 0, weight = 0;
			for (int t = 0; t < numTopics; t++) {
				weight = ((currentTypeTopicCounts.get(t) + corpus.gamma) / (corpus.citationsPerTopic[t] + corpus.gammaSum))
						* ((doc_topics.get(t) + corpus.alpha[t])/(doc.docLength+doc.citationSequence.getLength()-1+corpus.alphaSum));
				topicDistributionSum += weight;
				topicDistribution[t] = weight;
			}

			topic = random.nextDiscrete(topicDistribution, topicDistributionSum);

			doc.citationTopicAssignmentDuplicate[c]=topic;
			doc.citationTopicAssignment.put(keys[c], topic);   //no use now. 暂时不用此属性
			corpus.citationTopicCounts[keys[c]].adjustOrPutValue(topic, 1, 1);
			doc_topics.adjustOrPutValue(topic, 2, 1);
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

	/*public void estimateParameters(Corpus corpus) {
		// estimate parameters after every k iterations after burn-in
		numSamples++;
		// int[][] typeTopicCounts = new
		// int[corpus.vocabulary.size()][numTopics];
		// int[][] docTopicCounts = new int[corpus.docs.size()][numTopics];
		// int[][] citationTopicCounts = new int[corpus.docs.size()][numTopics];

		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			for (int i = 0; i < numTopics; i++) {
				if (numSamples > 1)
					corpus.theta_train[doc][i] *= (numSamples - 1);
				corpus.theta_train[doc][i] += (corpus.alpha[i] / (corpus.docs
						.get(doc).docLength + corpus.alphaSum));
			}
			int[] topics = corpus.docs.get(doc).topicAssignments.getFeatures();
			// System.out.println("=========");
			// System.out.println(topics.length);
			// System.out.println(corpus.docs.get(doc).topicAssignments.getLength());
			for (int i = 0; i < topics.length; i++) {
				corpus.theta_train[doc][topics[i]] += (1.0 / (corpus.docs
						.get(doc).docLength + corpus.alphaSum));
			}
			if (numSamples > 1) {
				for (int i = 0; i < numTopics; i++) {
					corpus.theta_train[doc][i] /= numSamples;
				}
			}
		}
		for (int k = 0; k < corpus.vocabulary.size(); k++) {
			for (int i = 0; i < numTopics; i++) {
				if (numSamples > 1)
					corpus.phi_train[k][i] *= (numSamples - 1);
				corpus.phi_train[k][i] += ((corpus.typeTopicCounts[k].get(i) + corpus.beta) / (corpus.tokensPerTopic[i] + corpus.betaSum));
				if (numSamples > 1)
					corpus.phi_train[k][i] /= numSamples;
			}
		}
		for (int k = 0; k < corpus.docAlphabet.size(); k++) {
			for (int i = 0; i < numTopics; i++) {
				if (numSamples > 1)
					corpus.psi_train[k][i] *= (numSamples - 1);
				corpus.psi_train[k][i] += ((corpus.citationTopicCounts[k]
						.get(i) + corpus.gamma) / (corpus.citationsPerTopic[i] + corpus.gammaSum));
				if (numSamples > 1)
					corpus.psi_train[k][i] /= (numSamples);
			}
		}

	}*/

	/**
	 * @author Echo
	 * @param corpus
	 * @return perplexity of the corpus ( not consider citations )
	 */
	public double modelPerplexity(Corpus corpus){
		double per = 0.0;
		double log_sum = 0.0;
		int word_sum = 0;
		for(int doc = 0; doc < corpus.docs.size(); doc++){
			ContextDocument cd = corpus.docs.get(doc);
			word_sum+=cd.docLength;
			TIntIntHashMap word_counts = corpus.docs.get(doc).word_counts;
			int[] words = word_counts.keys();
			//int[] words = corpus.docs.get(doc).wordSequence.getFeatures();
			//double sum_per_doc=0.0;
			for(int i=0; i<words.length; i++){
				int word_i = words[i];
				//System.out.print("word_i : "+word_i+"\t");
				double sum_per_word=0.0;
				for(int topic=0;topic<numTopics; topic++){
					sum_per_word += corpus.theta_train[doc][topic]*corpus.phi_train[word_i][topic]*word_counts.get(word_i);
					//if(corpus.theta_train[doc][topic]!=0&&corpus.phi_train[word_i][topic]!=0)
					//System.out.print("theta_train:"+corpus.theta_train[doc][topic]+", phi_train : "+corpus.phi_train[word_i][topic]);
				}
				//System.out.println("sum_per_word : "+sum_per_word);
				if(sum_per_word!=0)
					log_sum +=Math.log(sum_per_word);
				
				//sum_per_doc+=Math.log(sum_per_word);
			}
			//System.out.println("log_sum : "+log_sum);
			//log_sum = log_sum + Math.log(sum_per_doc);
		}
		//System.out.println("word_count_sum : "+word_sum);
		per = Math.exp(-log_sum/(word_sum+0.0));
		return per;
	}
	
	public double modelLogLikelihood(Corpus corpus) {
		double logLikelihood = 0.0;
		int nonZeroTopics;

		TIntIntHashMap[] typeTopicCounts = new TIntIntHashMap[corpus.numTypes];
		typeTopicCounts = corpus.typeTopicCounts;
		int[] tokensPerTopic = new int[numTopics];
		tokensPerTopic = corpus.tokensPerTopic;

		//int[] citationsPerTopic = new int[numTopics];
		//citationsPerTopic = corpus.citationsPerTopic;
		
		//out of index
		/*TIntIntHashMap[] citationTopicCounts = new TIntIntHashMap[corpus.numCitations+corpus.docAlphabet.size()];
		citationTopicCounts = corpus.citationTopicCounts;
		for (int i = 0; i < corpus.citationTopicCounts.length; i++)
			citationTopicCounts[i] = new TIntIntHashMap();
		for (int i = 0; i < corpus.typeTopicCounts.length; i++)
			typeTopicCounts[i] = new TIntIntHashMap();

		for (int doc_index = 0; doc_index < corpus.docs.size(); doc_index++) {
			ContextDocument doc = (ContextDocument) corpus.docs.get(doc_index);

			int[] topics = doc.topicAssignments.getFeatures();
			int[] words = doc.wordSequence.getFeatures();
			// put a check to see if wordlength exceeds words.length
			for (int i = 0; i < doc.sentenseLength; i++) {
				for (int j = doc.SentenseBoundaries.get(i).firstElement(); j < doc.SentenseBoundaries
						.get(i).lastElement(); j++) {

					if (doc.contextObject.contains(i)) {

						typeTopicCounts[words[j]].adjustOrPutValue(
								topics[j], 1, 1);
						tokensPerTopic[topics[j]]++;
						// for (int citations = 0; citations < doc.contextObject
						// .get(i).size(); citations++) {
						
														// alphabet
						// topics[j] = r.nextInt(numTopics);
						// corpus.typeTopicCounts[words[j]].adjustOrPutValue(
						// topics[j], 1, 1);
						citationTopicCounts[doc.contextObject.get(i)
								.get(0)].adjustOrPutValue(topics[j], 1, 1);
						// corpus.tokensPerTopic[topics[j]]++;
						citationsPerTopic[topics[j]]++;
						// }
					} else {

						typeTopicCounts[words[j]].adjustOrPutValue(
								topics[j], 1, 1);
						tokensPerTopic[topics[j]]++;
					}
				}
			}
		}
*/
		// The likelihood of the model is a combination of a
		// Dirichlet-multinomial for the words in each topic
		// and a Dirichlet-multinomial for the topics in each
		// document.

		// The likelihood function of a dirichlet multinomial is
		// Gamma( sum_i alpha_i ) prod_i Gamma( alpha_i + N_i )
		// prod_i Gamma( alpha_i ) Gamma( sum_i (alpha_i + N_i) )

		// So the log likelihood is
		// logGamma ( sum_i alpha_i ) - logGamma ( sum_i (alpha_i + N_i) ) +
		// sum_i [ logGamma( alpha_i + N_i) - logGamma( alpha_i ) ]

		// Do the documents first

		int[] topicCounts = new int[numTopics];
		double[] topicLogGammas = new double[numTopics];
		int[] docTopics;
		int sampleSize = 0;

		for (int topic = 0; topic < numTopics; topic++) {
			topicLogGammas[topic] = Dirichlet
					.logGammaStirling(corpus.alpha[topic]);
		}
		
		
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			LabelSequence topicSequence = (LabelSequence) corpus.docs.get(doc).topicAssignments;  //Never be initialized

			docTopics = topicSequence.getFeatures();
			sampleSize += docTopics.length;

			for (int token = 0; token < docTopics.length; token++) {
				topicCounts[docTopics[token]]++;
			}

			for (int topic = 0; topic < numTopics; topic++) {
				if (topicCounts[topic] > 0) {

					logLikelihood += (Dirichlet
							.logGammaStirling(corpus.alpha[topic]
									+ topicCounts[topic]) - topicLogGammas[topic]);
				}
			}

			// subtract the (count + parameter) sum term
			logLikelihood -= Dirichlet.logGammaStirling(corpus.alphaSum
					+ docTopics.length);

			Arrays.fill(topicCounts, 0);
		}

		// add the parameter sum term
		logLikelihood += corpus.docs.size()
				* Dirichlet.logGammaStirling(corpus.alphaSum);

		// And the topics

		// Count the number of type-topic pairs
		int nonZeroTypeTopics = 0;

		for (int type = 0; type < corpus.vocabulary.size(); type++) {
			int[] usedTopics = typeTopicCounts[type].keys();

			for (int topic : usedTopics) {
				int count = typeTopicCounts[type].get(topic);
				if (count > 0) {
					nonZeroTypeTopics++;
					logLikelihood += Dirichlet.logGammaStirling(corpus.beta
							+ count);
				}
			}
		}

		for (int topic = 0; topic < numTopics; topic++) {

			logLikelihood -= Dirichlet
					.logGammaStirling((corpus.beta * numTopics)
							+ tokensPerTopic[topic]);
		}

		logLikelihood += (Dirichlet.logGammaStirling(corpus.beta * numTopics))
				- (Dirichlet.logGammaStirling(corpus.beta) * nonZeroTypeTopics);

		return Math.exp(-1 * logLikelihood / sampleSize);
	}

	public double testPerplexity(Corpus corpus) {
		double ll = 0;
		int sampleSize = 0;
		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			ContextDocument document = (ContextDocument) corpus.docs.get(doc);
			sampleSize += document.docLength+document.citationTopicAssignmentDuplicate.length;
			TIntIntHashMap word_counts = corpus.docs.get(doc).word_counts;
			int[] keys = word_counts.keys();
			double tmp_ll = 0;
			// For words(content) part
			for (int i = 0; i < keys.length; i++) {
				for (int t = 0; t < numTopics; t++) {
					tmp_ll += corpus.theta_train[doc][t]
							* corpus.phi_train[keys[i]][t];
				}
				ll += (Math.log(tmp_ll) * word_counts.get(keys[i]));
			}
			// For citation part
			int[] citations = document.citationSequence.getFeatures();
			for (int j = 0; j < citations.length; j++){
				for(int t = 0; t < numTopics; t++){
					tmp_ll += corpus.theta_train[doc][t]
							* corpus.psi_train[citations[j]][t];
				}
				ll += Math.log(tmp_ll);
			}
		}
		//System.out.println("Samples=" + sampleSize);
		return Math.exp(-1 * ll / sampleSize);
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

	public static void main(String[] args) throws IOException{
		String input = "D:/WorkPaceEclipse_2/citation_LDA-master/dataset/VLDB_TS_Occur/paper_info.txt";
		String citationInput="D:/WorkPaceEclipse_2/citation_LDA-master/dataset/VLDB_TS_Occur/paper_citation.txt";
		int iterations = 5000;
		linkLDA_milestone_test linkLda = new linkLDA_milestone_test();
		System.out.println("Reading Data.....");
		Corpus corpus = linkLda.readData(input,citationInput);
		System.out.println("Done");
		//numTopics = 50;
		for(int i=0;i<=5;i++){
			linkLda.numTopics = 50;
			linkLda.InitializeParameters(corpus);
			LabelAlphabet labelAlp=linkLda.newLabelAlphabet(linkLda.numTopics);
		
		
		
			linkLda.InitializeAssignments(corpus, labelAlp);
			linkLda.sampleCorpus(corpus, iterations, true);
			linkLda.estimateParameters_single(corpus);
		
			System.out.println(iterations+" : "+linkLda.testPerplexity(corpus));
		}
	}
	
	public void estimateParameters_single(Corpus corpus) {
		for(int doc = 0; doc < corpus.docs.size();doc++){
			for(int i = 0;i < numTopics; i++) {
				corpus.theta_train[doc][i] = corpus.alpha[i]/(corpus.docs.get(doc).docLength+2*corpus.docs.get(doc).citationSequence.getLength()
						+corpus.alphaSum);
			}
			int[] topics = corpus.docs.get(doc).wordTopicAssignment;
			for(int i = 0;i < topics.length; i++){
				corpus.theta_train[doc][topics[i]] += 1.0/(corpus.docs.get(doc).docLength+2*corpus.docs.get(doc).citationSequence.getLength()
						+corpus.alphaSum);
			}
			topics = corpus.docs.get(doc).citationTopicAssignmentDuplicate;
			for(int i = 0;i < topics.length; i++){
				corpus.theta_train[doc][topics[i]] += 2.0/(corpus.docs.get(doc).docLength+2*corpus.docs.get(doc).citationSequence.getLength()
						+corpus.alphaSum);
			}
		}
		
		for (int k = 0; k < corpus.vocabulary.size(); k++) {
			for (int i = 0; i < numTopics; i++) {
				corpus.phi_train[k][i] = ((corpus.typeTopicCounts[k].get(i) + corpus.beta) / (corpus.tokensPerTopic[i] + corpus.betaSum));
			}
		}
		for (int k = 0; k < corpus.citationAlphabet.size(); k++) {
			for (int i = 0; i < numTopics; i++) {

				corpus.psi_train[k][i] = ((corpus.citationTopicCounts[k].get(i) + corpus.gamma) / (corpus.citationsPerTopic[i] + corpus.gammaSum));

			}
		}

	}
	
	/*public static void main(String[] args) throws FileNotFoundException, IOException {
		// TODO Auto-generated method stub
		// Read documents
		String input = "D:/WorkPaceEclipse_2/citation_LDA-master/dataset/VLDB_TS/paper_info.txt";
		String citationInput="D:/WorkPaceEclipse_2/citation_LDA-master/dataset/VLDB_TS/paper_citation.txt";
		String output = "D:/WorkPaceEclipse_2/citation_LDA-master/dataset/VLDB_TS/result_content_2citation_iter5000_beta_gamma_0.01_#2.txt";
		//String output = "D:/WorkPaceEclipse_2/citation_LDA-master/dataset/VLDB/test";
		
		int iterations = 5000;
		linkLDA_milestone_test linkLda = new linkLDA_milestone_test();
		System.out.println("Reading Data.....");
		Corpus corpus = linkLda.readData(input,citationInput);
		System.out.println("Done");
		//numTopics = 50;
		linkLda.numTopics = 50;
		linkLda.InitializeParameters(corpus);
	
		LabelAlphabet labelAlp=linkLda.newLabelAlphabet(linkLda.numTopics);
		linkLda.InitializeAssignments(corpus,labelAlp);
		linkLda.sampleCorpus(corpus, iterations, true);
		linkLda.estimateParameters_single(corpus);
		
		
		for(int i=0;i<linkLda.numTopics;i++)
		{
			System.out.println(corpus.alpha[i]);
		}
		
	
		for(int j=0;j<corpus.docs.size();j++)
		{
			int[] topics = corpus.docs.get(j).topicAssignments.getFeatures();
			System.out.println("Doc"+corpus.docAlphabet.lookupObject(j));
			for (int k = 0; k <corpus.docs.get(j).wordSequence.getLength(); k++)
			{
				System.out.print(corpus.vocabulary.lookupObject(corpus.docs.get(j).wordSequence.getFeatures()[k])+" "+topics[k]);
				
			}
			
			System.out.println();
		}
		
		ArrayList<TreeSet<IDSorter>> topicSortedWords = linkLda.getSortedWords(corpus,linkLda.numTopics);
		System.out.println("================================");
		Formatter out = new Formatter(new StringBuilder(), Locale.US);
		BufferedWriter bw1= new BufferedWriter(new OutputStreamWriter(new FileOutputStream(new File(output),true)));
		bw1.write("# Topic_word");  
		bw1.newLine();
		for (int topic = 0; topic < linkLda.numTopics; topic++) {
			Iterator<IDSorter> iterator = topicSortedWords.get(topic).iterator();
			
			out = new Formatter(new StringBuilder(), Locale.US);
			out.format ("%d\t",topic);
			int rank = 0;
			while (iterator.hasNext() && rank<20) {
				IDSorter idCountPair = iterator.next();
				out.format("%s (%.4f) ", corpus.vocabulary.lookupObject(idCountPair.getID()), idCountPair.getWeight());
				rank++;
			}
		  System.out.println(out);
			String line=out.toString();
			bw1.write(line);  
			bw1.newLine();
		}
		
		bw1.write("# Topic_Citation");  
		bw1.newLine();
		HashMap<Integer,String> citationMap=new HashMap<Integer,String>();
		ArrayList<TreeSet<IDSorter>> topicSortedCitation = linkLda.getSortedCitation(corpus,linkLda.numTopics,citationMap);
		double sum=0;
	
		for (int topic = 0; topic < linkLda.numTopics; topic++) {
			Iterator<IDSorter> iterator = topicSortedCitation.get(topic).iterator();
			
			out = new Formatter(new StringBuilder(), Locale.US);
			out.format ("%d\n",topic);
			int rank = 0;
			while (iterator.hasNext()) {
				IDSorter idCountPair = iterator.next();
				if(idCountPair.getWeight()>0.0001)
				{
				String title=null;
				int year=-1;
				for(int i=0;i<fullpaperList.size();i++)
				{
					if(fullpaperList.get(i).getID().equals(citationMap.get(idCountPair.getID())))
					{
						title=fullpaperList.get(i).getTitle();
						year = fullpaperList.get(i).getYear();
					break;
					}
				}
				out.format("%s\t%.4f\t%s\t%d\n", citationMap.get(idCountPair.getID()), idCountPair.getWeight(),title,year);
				rank++;
			}
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
					year = fullpaperList.get(z).getYear();
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

*/	
}
