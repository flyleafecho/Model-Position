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
import java.util.Formatter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Locale;
import java.util.TreeSet;

import edu.psu.types.ContextDocument;
import edu.psu.types.Corpus;
import edu.psu.types.IDSorter;
import edu.psu.types.LabelAlphabet;
import edu.psu.types.Paper;
import edu.psu.util.Randoms;
import gnu.trove.TIntIntHashMap;

public class CitePart_LDA extends Model {
	public int numParts = 0;
	public ArrayList<Paper> paperList = new ArrayList<Paper>();
	public ArrayList<Paper> fullpaperList = new ArrayList<Paper>();

	public Corpus readData(String docDir, String citationDir) throws IOException, FileNotFoundException{			
		File paper_info = new File(docDir);
		BufferedReader metaReader = new BufferedReader(new InputStreamReader(new FileInputStream(paper_info),"utf-8"));
		File paper_citation = new File(citationDir);
		BufferedReader citationReader = new BufferedReader(new InputStreamReader(new FileInputStream(paper_citation),"utf-8"));
		
		ArrayList<String> papers = new ArrayList<String>();
		ArrayList<String> citations = new ArrayList<String>();
		
		String lineMeta = null;
		String lineCitation = null;
		
		int lineNum = 0;
		while((lineCitation=citationReader.readLine())!=null){
			lineNum++;
			if(lineNum%2==1){
				papers.add(lineCitation);
			}else{
				citations.add(lineCitation);
			}
		}
		citationReader.close();
		
		while((lineMeta=metaReader.readLine())!=null){
			String[] parts = lineMeta.split("\\t");
			Paper paper = new Paper();
			paper.setID(parts[0]);
			paper.setTitle(parts[1]);
			paper.setYear(Integer.parseInt(parts[2]));
			fullpaperList.add(paper);
			
			if(papers.contains(parts[0])){
				paper.setAbs(parts[3]);
				paperList.add(paper);
			}
		}
		metaReader.close();
		
		Corpus corpus = new Corpus(paperList.size());
		for (int i = 0; i < paperList.size(); i++) {
			if (!Corpus.docAlphabet.contains(paperList.get(i).getID())) {
				Corpus.docAlphabet.lookupIndex(paperList.get(i).getID());
			}
		}
		
		for (int i = 0; i < paperList.size(); i++) {
			ContextDocument doc = new ContextDocument();
			doc.add_milestone(paperList.get(i).getID(), Corpus.docAlphabet.lookupIndex(paperList.get(i).getID()),(paperList.get(i).getTitle()+" "+paperList.get(i).getAbs()).toLowerCase(),papers,citations);
			corpus.addDocument(doc,Corpus.docAlphabet.lookupIndex(paperList.get(i).getID()));
		}
		
		System.out.println("Total Cited Documents:"+ Corpus.citationAlphabet.size());
		System.out.println("Total Documents:" + corpus.docs.size());
		System.out.println("Total Vcabulary Size:"+ Corpus.vocabulary.size());
		return corpus;
	}
	
	public void InitializeParameters(Corpus corpus){
		corpus.numCitations = Corpus.citationAlphabet.size();
		corpus.numTypes = Corpus.vocabulary.size();
		corpus.numParts = 2;
		numParts = 2;
		
		this.alpha = new double[numTopics];
		this.alphaSum = 0;
		for (int i = 0; i < numTopics; i++) {
			this.alpha[i] = 50.0 / (double) numTopics;
			this.alphaSum += this.alpha[i];
		}
		//beta is the prior distribution for topic_citation
		this.beta = 0.01;
		this.betaSum = this.beta * corpus.numCitations;
		//gamma is the prior distribution for doc_part
		this.gamma = 0.01;
		this.gammaSum = this.gamma * corpus.numParts;
		
		//corpus.typeTopicCounts = new TIntIntHashMap[corpus.numTypes];
		//corpus.tokensPerTopic = new int[numTopics];
		corpus.citationsPerTopic = new int[numTopics];
		corpus.citationsPerTopicByPart = new int[numParts][numTopics];
		
		this.phi_train = new double[Corpus.citationAlphabet.size()][numTopics];	//topic_citation
		this.theta_train = new double[Corpus.docAlphabet.size()][numTopics];	//doc_topic
		
		this.theta_train_part = new double[numParts][Corpus.docAlphabet.size()][numTopics];
		
		this.psi_train = new double[Corpus.docAlphabet.size()][numParts];	//doc_part
		//citationTopicCounts(大小不加docAlphabet应该也是对的)
		//Array, per obj is a map which represents specific citation's times are assigned to topics
		corpus.citationTopicCounts = new TIntIntHashMap[corpus.numCitations];
		for (int i = 0; i < corpus.citationTopicCounts.length; i++)
			corpus.citationTopicCounts[i] = new TIntIntHashMap();
		/*for (int i = 0; i < corpus.typeTopicCounts.length; i++)
			corpus.typeTopicCounts[i] = new TIntIntHashMap();*/
	}
	
	public void InitializeAssignments(Corpus corpus, LabelAlphabet topicAlphabet){
		this.random = new Randoms();
		for(int doc_index = 0; doc_index < corpus.docs.size(); doc_index++){
			ContextDocument doc = corpus.docs.get(doc_index);
			//doc.citationTopicAssignmentDuplicate = new int[doc.citationSequence.getLength()];
			doc.citationTAByPart = new int[numParts][];
			for(int i = 0; i < numParts; i++){
				doc.citationTAByPart[i] = new int[doc.citationSequenceByPart[i].getLength()];
			}
			
			//Initiate different parts separately
			for(int i = 0; i < numParts; i++){
				int[] keys = doc.citationSequenceByPart[i].getFeatures();
				for(int j = 0; j < doc.citationSequenceByPart[i].getLength(); j++){
					int topic = r.nextInt(numTopics);
					doc.citationTAByPart[i][j] = topic;
					
					corpus.citationTopicCounts[keys[j]].adjustOrPutValue(topic, 1, 1);
					corpus.citationsPerTopicByPart[i][topic]++;			//used for estimate doc_topic
				}
			}
			corpus.docs.set(doc_index, doc);
		}
	}
	
	public void sampleCorpus(Corpus corpus, int iterations, boolean isContextAware){
		for(int i = 0; i < iterations; i++){
			for(int j = 0; j < corpus.docs.size(); j++){
				sampleOneDocument(corpus,(ContextDocument) corpus.docs.get(j));
			}
		}
	}
	
	public void sampleOneDocument(Corpus corpus, ContextDocument doc){
		
		TIntIntHashMap[] doc_topic_by_part = new TIntIntHashMap[numParts];  //topic, and the occurrance time of the topic
		for(int i = 0; i < numParts; i++)
			doc_topic_by_part[i] = new TIntIntHashMap();
		
		for(int i = 0; i < numParts; i++){
			//int[] keys = doc.citationSequenceByPart[i].getFeatures();
			for (int t = 0; t < doc.citationSequenceByPart[i].getLength(); t++) 
			{
				doc_topic_by_part[i].adjustOrPutValue(doc.citationTAByPart[i][t], 1, 1);
			}
		}

		for(int i = 0; i < numParts; i++){
			int[] keys = doc.citationSequenceByPart[i].getFeatures();	//keys' size are bigger than or equal to the length of FS 
			for(int t = 0; t < doc.citationSequenceByPart[i].getLength(); t++){
				int topic = doc.citationTAByPart[i][t];
				
				corpus.citationsPerTopicByPart[i][topic]--;
				corpus.citationTopicCounts[keys[t]].adjustOrPutValue(topic, -1, 0);
				doc_topic_by_part[i].adjustOrPutValue(topic, -1, 0);
				
				//TIntIntHashMap currentTypeTopicCounts = corpus.citationTopicCounts[keys[t]];
				double[] topicDistribution = new double[numTopics];
				double topicDistributionSum = 0, weight = 0;
								
				for(int k = 0; k < numTopics; k++){
					int topicCitationSum = corpus.citationsPerTopicByPart[0][k] + corpus.citationsPerTopicByPart[1][k];
					weight = (doc.citationSequenceByPart[i].getLength()+gamma - 1)*(doc_topic_by_part[i].get(k) + alpha[k])*
								(corpus.citationTopicCounts[keys[t]].get(k)+beta)/(topicCitationSum + betaSum);
					topicDistributionSum += weight;
					topicDistribution[k] = weight;
				}
				
				topic = random.nextDiscrete(topicDistribution,topicDistributionSum);
				doc.citationTAByPart[i][t] = topic;
				
				corpus.citationsPerTopicByPart[i][topic]++;
				corpus.citationTopicCounts[keys[t]].adjustOrPutValue(topic, 1, 1);
				doc_topic_by_part[i].adjustOrPutValue(topic, 1, 1);
			}
		}
	}
	
	public void estimateParameters_single(Corpus corpus){
		for(int i = 0; i < numParts; i++){
			for(int d = 0; d < corpus.docs.size(); d++){
				this.psi_train[d][i] = (corpus.docs.get(d).citationSequenceByPart[i].getLength() + gamma)/
											(corpus.docs.get(d).citationSequence.getLength() + gammaSum);	//this is the final value of psi
				for(int k = 0; k < numTopics; k++){
					this.theta_train_part[i][d][k] = alpha[k]/(corpus.docs.get(d).citationSequenceByPart[i].getLength() + alphaSum);
				}
				
				ContextDocument doc = corpus.docs.get(d);
				int[] topics = doc.citationTAByPart[i];
				for(int k : topics){
					this.theta_train_part[i][d][k] += 1/(corpus.docs.get(d).citationSequenceByPart[i].getLength() + alphaSum); 
				}
			}
		}
		
		for(int c = 0; c < corpus.citationAlphabet.size(); c++){
			for(int k = 0; k < numTopics; k++){
				int topicCitationSum = corpus.citationsPerTopicByPart[0][k] + corpus.citationsPerTopicByPart[1][k];
				this.phi_train[c][k] = (corpus.citationTopicCounts[c].get(k) + beta)/(topicCitationSum+betaSum);
			}
		}
	}
	
	public void writeResult(String outDir, Corpus corpus) throws IOException{
		File output = new File(outDir);
		Formatter out = new Formatter(new StringBuilder(), Locale.US);
		BufferedWriter bw= new BufferedWriter(new OutputStreamWriter(new FileOutputStream(output,false)));
		
		bw.write("# Topic_Citation");  
		bw.newLine();
		
		HashMap<Integer,String> citationMap=new HashMap<Integer,String>();
		ArrayList<TreeSet<IDSorter>> topicSortedCitation = getSortedCitation(corpus,numTopics,citationMap);
		
		for(int topic = 0; topic < numTopics; topic++){
			Iterator<IDSorter> iterator = topicSortedCitation.get(topic).iterator();
			
			out = new Formatter(new StringBuilder(), Locale.US);
			out.format ("%d\n",topic);
			int rank = 0;
			while (iterator.hasNext()) {
				IDSorter idCountPair = iterator.next();
				if(idCountPair.getWeight()>0.0001){
					String title=null;
					int year=-1;
					for(int i = 0; i < fullpaperList.size(); i++){
						if(fullpaperList.get(i).getID().equals(citationMap.get(idCountPair.getID()))){
							title=fullpaperList.get(i).getTitle();
							year = fullpaperList.get(i).getYear();
							break;
						}
					}
					out.format("%s\t%.4f\t%s\t%d\n", citationMap.get(idCountPair.getID()), idCountPair.getWeight(),title,year);
					rank++;
				}
			}
			String line = out.toString();
			System.out.println(line);
			bw.write(line);
		}
		
		bw.write("# Doc_Topic");  
		bw.newLine();
		for(int i = 0; i < corpus.docAlphabet.size(); i++){
			out = new Formatter(new StringBuilder(), Locale.US);
			String id = (String) corpus.docAlphabet.lookupObject(i);
			String title=null;
			int year=-1;
			for(int z = 0; z < fullpaperList.size(); z++){
				if(fullpaperList.get(z).getID().equals(id)){
					title=fullpaperList.get(z).getTitle();
					year = fullpaperList.get(z).getYear();
					break;
				}
			}
			out.format("%s\t", corpus.docAlphabet.lookupObject(i));
			double[] max ={0,0}; int[] assign = {-1,-1};
			for(int part = 0; part < numParts; part++){
				for(int topic = 0; topic < numTopics; topic++){
					if(theta_train_part[part][i][topic]>max[part]){
						max[part] = theta_train_part[part][i][topic];
						assign[part]=topic;
					}
					out.format ("%.4f\t",theta_train_part[part][i][topic]);
				}
			}
			out.format("%d\t%.4f\t%d\t%.4f\t%s\t%d", assign[0],max[0],assign[1],max[1],title,year);
			String line = out.toString();
			System.out.println(line);
			bw.write(line);
			bw.newLine();
		}
		
		out.close();
		bw.close();
	}
	
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
				topicSortedCitation.get(i).add(new IDSorter(index, phi_train[type][i]));
			}
		}
		return topicSortedCitation;
	}
	
	public double sampleLikelihood(int numSamples, Corpus corus){return 0;}
	
	public double sampleLikelihood(int numSamples, Corpus corpus,
			TIntIntHashMap split){return 0;}
	
	public static void main(String[] args) throws Exception{
		String paper_info = "D:\\WorkPaceEclipse_2\\citation_LDA-master\\dataset\\VLDB_TS_CP\\paper_info.txt";
		String paper_citation = "D:\\WorkPaceEclipse_2\\citation_LDA-master\\dataset\\VLDB_TS_CP\\paper_citation_part.txt";
		String result = "D:\\WorkPaceEclipse_2\\citation_LDA-master\\dataset\\VLDB_TS_CP\\out.txt";
		
		int iterations = 5000;
		
		CitePart_LDA cpLDA = new CitePart_LDA();
		
		System.out.println("Reading Data......");
		Corpus corpus = cpLDA.readData(paper_info, paper_citation);
		
		cpLDA.numTopics = 50;
		cpLDA.InitializeParameters(corpus);
		
		LabelAlphabet labelAlp= null;	//Just for matching the abstract method in Model
		cpLDA.InitializeAssignments(corpus,labelAlp);
		
		cpLDA.sampleCorpus(corpus, iterations, false);
		cpLDA.estimateParameters_single(corpus);
		
		cpLDA.writeResult(result,corpus);
	}
}
