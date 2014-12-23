package edu.psu.topic;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;
import java.util.TreeSet;

import edu.psu.types.*;
import edu.psu.util.Randoms;
import gnu.trove.TIntIntHashMap;

public abstract class Model {
	public abstract void InitializeParameters(Corpus corpus);

	public abstract void InitializeAssignments(Corpus corpus,
			LabelAlphabet topicAlphabet);

	public  int numTopics;
	//public int numSamples = 0;
	public Random r = new Random();
	public Randoms random;
	
	//Model Parameters
	public double[] alpha;	//Prior distribution of doc_topic
	public double alphaSum;
	public double beta;		//Prior distribution of topic_citation
	public double betaSum;
	
	public double gamma;	//Prior distribution of doc_part-weight
	public double gammaSum;
	
	//Model Results
	public double[][] theta_train;	//M*K
	public double[][] phi_train;	//Vocabulary*K
	public double[][] psi_train;	//M*2
	
	// for dirichlet estimation
	public int[] docLengthCounts; // histogram of document sizes
	public int[][] topicDocCounts;

	public abstract void sampleOneDocument(Corpus corpus, ContextDocument doc);

	public abstract void sampleCorpus(Corpus corpus, int iterations,
			boolean isContextAware);
	
	public void learnParameters(Corpus corpus){
		docLengthCounts = new int[corpus.maxTokens+1];
		topicDocCounts = new int[numTopics][corpus.maxTokens + 1];
		for (int d = 0; d < corpus.docs.size(); d++) {
			ContextDocument doc = (ContextDocument) corpus.docs.get(d);
			
			docLengthCounts[doc.wordLength]++;
			int[] topics = doc.topicAssignments.getFeatures();
			int[] topicCounts=  new int[numTopics];
			for(int i=0;i<topics.length;i++){
				topicCounts[topics[i]]++;
				
			}
			for (int topic=0; topic < numTopics; topic++) {
				
				topicDocCounts[topic][ topicCounts[topic] ]++;
				//System.out.print(topicCounts[topic] + "\t");
			}
			//System.out.println();
		}
		alphaSum = Dirichlet.learnParameters(alpha, topicDocCounts, docLengthCounts);
	}

	
	public abstract void estimateParameters_single(Corpus corpus);
	/**
	 *  Return an array of sorted sets (one set per topic). Each set 
	 *   contains IDSorter objects with integer keys into the alphabet.
	 *   To get direct access to the Strings, use getTopWords().
	 */
	public ArrayList<TreeSet<IDSorter>> getSortedWords (Corpus corpus,int numTopics) {
	
		ArrayList<TreeSet<IDSorter>> topicSortedWords = new ArrayList<TreeSet<IDSorter>>(numTopics);

		// Initialize the tree sets
		for (int topic = 0; topic < numTopics; topic++) {
			topicSortedWords.add(new TreeSet<IDSorter>());
		}

		// Collect counts
		for (int type = 0; type < corpus.vocabulary.size(); type++) {
			for (int i = 0; i < numTopics; i++) {
				phi_train[type][i] = ((corpus.typeTopicCounts[type].get(i) + beta) / (corpus.tokensPerTopic[i] + betaSum));
				topicSortedWords.get(i).add(new IDSorter(type, phi_train[type][i]));
			}

		}

		return topicSortedWords;
	}
	
	// Need overwrited in specific Models
	public ArrayList<TreeSet<IDSorter>> getSortedCitation (Corpus corpus,int numTopics,HashMap<Integer,String> citationMap) {
		
		ArrayList<TreeSet<IDSorter>> topicSortedCitation = new ArrayList<TreeSet<IDSorter>>(numTopics);

		// Initialize the tree sets
		for (int topic = 0; topic < numTopics; topic++) {
			topicSortedCitation.add(new TreeSet<IDSorter>());
		}

		// Collect counts
		ArrayList<String> entries=Corpus.citationAlphabet.entries;
		for (int type = 0; type < Corpus.citationAlphabet.size(); type++) {
			int index=Corpus.citationAlphabet.lookupIndex(entries.get(type));
			citationMap.put(index, entries.get(type));
			for (int i = 0; i < numTopics; i++) {			
				psi_train[type][i] = ((corpus.citationTopicCounts[index].get(i) + gamma) / (corpus.citationsPerTopic[i] + gammaSum));
				topicSortedCitation.get(i).add(new IDSorter(index, psi_train[type][i]));
			}

		}

		return topicSortedCitation;
	}
	
	public abstract double sampleLikelihood(int numSamples, Corpus corus);

	public abstract double sampleLikelihood(int numSamples, Corpus corpus,
			TIntIntHashMap split);

	public double testPerplexity(int numSamples, Corpus corpus) {
		double ll = 0;
		int sampleSize = 0, ss = 0;

		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			ContextDocument document = (ContextDocument) corpus.docs.get(doc);
			sampleSize += document.wordLength;
			TIntIntHashMap word_counts = corpus.docs.get(doc).word_counts;
			int[] keys = word_counts.keys();
			double tmp_ll = 0;
			for (int i = 0; i < keys.length; i++) {
				for (int t = 0; t < numTopics; t++) {
					tmp_ll += theta_train[doc][t]
							* phi_train[keys[i]][t];

				}
				ll += (Math.log(tmp_ll) * word_counts.get(keys[i]));
				ss += word_counts.get(keys[i]);

			}

		}
		System.out.println("Samples=" + sampleSize + "docs:"
				+ corpus.docs.size() + "ss:" + ss);

		// return Math.exp(-1 * ll / sampleSize);
		return ll;
	}

	public double modelLogLikelihood(Corpus corpus) {
		double logLikelihood = 0.0;
		int nonZeroTopics;

		TIntIntHashMap[] typeTopicCounts = new TIntIntHashMap[corpus.numTypes];
		int[] tokensPerTopic = new int[numTopics];
		int[] citationsPerTopic = new int[numTopics];

		TIntIntHashMap[] citationTopicCounts = new TIntIntHashMap[corpus.numCitations];
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

						typeTopicCounts[words[j]].adjustOrPutValue(topics[j],
								1, 1);
						tokensPerTopic[topics[j]]++;
						// for (int citations = 0; citations < doc.contextObject
						// .get(i).size(); citations++) {

						// alphabet
						// topics[j] = r.nextInt(numTopics);
						// corpus.typeTopicCounts[words[j]].adjustOrPutValue(
						// topics[j], 1, 1);
						citationTopicCounts[doc.contextObject.get(i).get(0)]
								.adjustOrPutValue(topics[j], 1, 1);
						// corpus.tokensPerTopic[topics[j]]++;
						citationsPerTopic[topics[j]]++;
						// }
					} else {

						typeTopicCounts[words[j]].adjustOrPutValue(topics[j],
								1, 1);
						tokensPerTopic[topics[j]]++;
					}
				}
			}
		}

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
					.logGammaStirling(alpha[topic]);
		}

		for (int doc = 0; doc < corpus.docs.size(); doc++) {
			LabelSequence topicSequence = (LabelSequence) corpus.docs.get(doc).topicAssignments;

			docTopics = topicSequence.getFeatures();
			sampleSize += docTopics.length;

			for (int token = 0; token < docTopics.length; token++) {
				topicCounts[docTopics[token]]++;
			}

			for (int topic = 0; topic < numTopics; topic++) {
				if (topicCounts[topic] > 0) {

					logLikelihood += (Dirichlet
							.logGammaStirling(alpha[topic]
									+ topicCounts[topic]) - topicLogGammas[topic]);
				}
			}

			// subtract the (count + parameter) sum term
			logLikelihood -= Dirichlet.logGammaStirling(alphaSum
					+ docTopics.length);

			Arrays.fill(topicCounts, 0);
		}

		// add the parameter sum term
		logLikelihood += corpus.docs.size()
				* Dirichlet.logGammaStirling(alphaSum);

		// And the topics

		// Count the number of type-topic pairs
		int nonZeroTypeTopics = 0;

		for (int type = 0; type < corpus.vocabulary.size(); type++) {
			int[] usedTopics = typeTopicCounts[type].keys();

			for (int topic : usedTopics) {
				int count = typeTopicCounts[type].get(topic);
				if (count > 0) {
					nonZeroTypeTopics++;
					logLikelihood += Dirichlet.logGammaStirling(beta
							+ count);
				}
			}
		}

		for (int topic = 0; topic < numTopics; topic++) {

			logLikelihood -= Dirichlet
					.logGammaStirling((beta * numTopics)
							+ tokensPerTopic[topic]);
		}

		logLikelihood += (Dirichlet.logGammaStirling(beta * numTopics))
				- (Dirichlet.logGammaStirling(beta) * nonZeroTypeTopics);

		//return Math.exp(-1 * logLikelihood / sampleSize);
		return logLikelihood;
	}

	public static Model factory(String type) {

		if (type.toLowerCase().equals("citelda")) {
			return new citeLDA();
		} else if (type.toLowerCase().equals("linklda")) {
			return new linkLDA();
		} else if (type.toLowerCase().equals("citeplsalda")) {
			return new citePlsaLDA();
		} else if (type.toLowerCase().equals("linkplsalda")) {
			return new linkPlsaLDA();
		} else
			return null;
	}

}
