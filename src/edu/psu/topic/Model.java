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
	public double[][][] theta_train_part;	//P*M*K
	public double[][] phi_train;	//Vocabulary*K
	public double[][] psi_train;	//M*2
	
	// for dirichlet estimation
	public int[] docLengthCounts; // histogram of document sizes
	public int[][] topicDocCounts;

	public abstract void sampleOneDocument(Corpus corpus, ContextDocument doc);

	public abstract void sampleCorpus(Corpus corpus, int iterations,
			boolean isContextAware);
	
	public abstract void estimateParameters_single(Corpus corpus);
	/**
	 *  Return an array of sorted sets (one set per topic). Each set 
	 *   contains IDSorter objects with integer keys into the alphabet.
	 *   To get direct access to the Strings, use getTopWords().
	 */

	public abstract double sampleLikelihood(int numSamples, Corpus corus);

	public abstract double sampleLikelihood(int numSamples, Corpus corpus,
			TIntIntHashMap split);

}
