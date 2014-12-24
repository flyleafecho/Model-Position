package edu.psu.topic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

import edu.psu.types.ContextDocument;
import edu.psu.types.Corpus;
import edu.psu.types.LabelAlphabet;
import edu.psu.types.Paper;
import edu.psu.util.Randoms;
import gnu.trove.TIntIntHashMap;

public class CitePart_LDA extends Model {
	public int numParts = 0;

	public Corpus readData(String docDir, String citationDir) throws IOException, FileNotFoundException{
		ArrayList<Paper> paperList = new ArrayList<Paper>();
		ArrayList<Paper> fullpaperList = new ArrayList<Paper>();
		
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
		this.phi_train = new double[Corpus.citationAlphabet.size()][numTopics];	//topic_citation
		this.theta_train = new double[Corpus.docAlphabet.size()][numTopics];	//doc_topic
		this.psi_train = new double[Corpus.docAlphabet.size()][numParts];	//doc_part
		//citationTopicCounts(大小不加docAlphabet应该也是对的)
		//Array, per obj is a map which represents specific citation's times are assigned to topics
		corpus.citationTopicCounts = new TIntIntHashMap[corpus.numCitations];
		for (int i = 0; i < corpus.citationTopicCounts.length; i++)
			corpus.citationTopicCounts[i] = new TIntIntHashMap();
		/*for (int i = 0; i < corpus.typeTopicCounts.length; i++)
			corpus.typeTopicCounts[i] = new TIntIntHashMap();*/
	}
	
	public void InitializeAssignments(Corpus corpus,
			LabelAlphabet topicAlphabet){
		this.random = new Randoms();
		for(int doc_index = 0; doc_index < corpus.docs.size(); doc_index++){
			ContextDocument doc = corpus.docs.get(doc_index);
		}
	}
	public static void main(String[] args) throws Exception{
		String paper_info = "";
		String paper_citation = "";
		
		int iteration = 5000;
		
		CitePart_LDA cpLDA = new CitePart_LDA();
		
		System.out.println("Reading Data......");
		Corpus corpus = cpLDA.readData(paper_info, paper_citation);
		
		cpLDA.numTopics = 50;
		cpLDA.InitializeParameters(corpus);
		
		LabelAlphabet labelAlp= null;	//Just for matching the abstract method in Model
		cpLDA.InitializeAssignments(corpus,labelAlp);
	}
}
