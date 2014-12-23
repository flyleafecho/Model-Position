package test;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;

public class Test {
	public static void main(String[] args) throws Exception{
		File file = new File("D:/WorkPaceEclipse_2/citation_LDA-master/dataset/test/id.txt");
		File fout = new File("D:/WorkPaceEclipse_2/citation_LDA-master/dataset/test/out.txt");
		BufferedReader br = new BufferedReader(new FileReader(file));
		FileWriter fw = new FileWriter(fout);
		String str = null;
		while((str=br.readLine())!=null){
			String out = str.substring(0,str.indexOf("(")).trim();
			fw.write(out+"\n");
		}
		br.close();
		fw.close();
	}
}
