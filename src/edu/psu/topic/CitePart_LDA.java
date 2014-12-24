package edu.psu.topic;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;

import edu.psu.types.Paper;

public class CitePart_LDA extends Model {
	public Corpus readData(String docDir, String citationDir) throws IOException, FileNotFoundException{
		ArrayList<Paper> paperList = new ArrayList<Paper>();
		ArrayList<Paper> fullpaperList = new ArrayList<Paper>();
		
		File paper_info = new File(docDir);
		BufferedReader metareader = new BufferedReader(new InputStreamReader(new FileInputStream(paper_info),"utf-8"));
	}
	
	public static void main(String[] args){
		
	}
}
