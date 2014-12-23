package edu.psu.util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;


public class RemoveStopWord {
	HashSet<String> stoplist = null;
	boolean caseSensitive = true;
   String [] tokens;
	private HashSet<String> newDefaultStopList ()
	{
		HashSet<String> sl = new HashSet<String>();
		for (int i = 0; i < stopwords.length; i++)
			sl.add (stopwords[i]);
		return sl;
	}
	public RemoveStopWord(File stoplistFile, String encoding, String[] tokens,boolean includeDefault,boolean caseSensitive) {
       if (! includeDefault) { stoplist = new HashSet<String>(); }
        else { stoplist = newDefaultStopList(); }

           addStopWords (fileToStringArray(stoplistFile, encoding));
           this.tokens=tokens;
             this.caseSensitive = caseSensitive;
}
	  public void addStopWords (String[] words)
	   {
		for (int i = 0; i < words.length; i++)
			stoplist.add (words[i]);
	   }

	
	
	private String[] fileToStringArray (File f, String encoding)
	{
		ArrayList<String> wordarray = new ArrayList<String>();

		try {

			BufferedReader input = null;
			if (encoding == null) {
				input = new BufferedReader (new FileReader (f));
			}
			else {
				input = new BufferedReader( new InputStreamReader( new FileInputStream(f), encoding ));
			}
			String line;

			while (( line = input.readLine()) != null) {
				String[] words = line.split ("\\s+");
				for (int i = 0; i < words.length; i++)
					wordarray.add (words[i]);
			}

		} catch (IOException e) {
			throw new IllegalArgumentException("Trouble reading file "+f);
		}
		return (String[]) wordarray.toArray(new String[]{});
	}
	
	public String [] RemoveStopWords()
	{
		ArrayList<String> array=new ArrayList<String>();
		String [] tokenRemove;
		for(int i=0;i<tokens.length;i++)
		{
			if (! stoplist.contains (caseSensitive ? tokens[i] : tokens[i].toLowerCase())) {
				// xxx Should we instead make and add a copy of the Token?
				array.add(tokens[i]);
			}	 
		}
		tokenRemove=(String[]) array.toArray(new String[]{});
		return tokenRemove;
	}
	static final String[] stopwords =
		{
			"a",
			"able",
			"about",
			"above",
			"according",
			"accordingly",
			"across",
			"actually",
			"after",
			"afterwards",
			"again",
			"against",
			"all",
			"allow",
			"allows",
			"almost",
			"alone",
			"along",
			"already",
			"also",
			"although",
			"always",
			"am",
			"among",
			"amongst",
			"an",
			"and",
			"another",
			"any",
			"anybody",
			"anyhow",
			"anyone",
			"anything",
			"anyway",
			"anyways",
			"anywhere",
			"apart",
			"appear",
			"appreciate",
			"appropriate",
			"are",
			"around",
			"as",
			"aside",
			"ask",
			"asking",
			"associated",
			"at",
			"available",
			"away",
			"awfully",
			"b",
			"be",
			"became",
			"because",
			"become",
			"becomes",
			"becoming",
			"been",
			"before",
			"beforehand",
			"behind",
			"being",
			"believe",
			"below",
			"beside",
			"besides",
			"best",
			"better",
			"between",
			"beyond",
			"both",
			"brief",
			"but",
			"by",
			"c",
			"came",
			"can",
			"cannot",
			"cant",
			"cause",
			"causes",
			"certain",
			"certainly",
			"changes",
			"clearly",
			"co",
			"com",
			"come",
			"comes",
			"concerning",
			"consequently",
			"consider",
			"considering",
			"contain",
			"containing",
			"contains",
			"corresponding",
			"could",
			"course",
			"currently",
			"d",
			"definitely",
			"described",
			"despite",
			"did",
			"different",
			"do",
			"does",
			"doing",
			"done",
			"down",
			"downwards",
			"during",
			"e",
			"each",
			"edu",
			"eg",
			"eight",
			"either",
			"else",
			"elsewhere",
			"enough",
			"entirely",
			"especially",
			"et",
			"etc",
			"even",
			"ever",
			"every",
			"everybody",
			"everyone",
			"everything",
			"everywhere",
			"ex",
			"exactly",
			"example",
			"except",
			"f",
			"far",
			"few",
			"fifth",
			"first",
			"five",
			"followed",
			"following",
			"follows",
			"for",
			"former",
			"formerly",
			"forth",
			"four",
			"from",
			"further",
			"furthermore",
			"g",
			"get",
			"gets",
			"getting",
			"given",
			"gives",
			"go",
			"goes",
			"going",
			"gone",
			"got",
			"gotten",
			"greetings",
			"h",
			"had",
			"happens",
			"hardly",
			"has",
			"have",
			"having",
			"he",
			"hello",
			"help",
			"hence",
			"her",
			"here",
			"hereafter",
			"hereby",
			"herein",
			"hereupon",
			"hers",
			"herself",
			"hi",
			"him",
			"himself",
			"his",
			"hither",
			"hopefully",
			"how",
			"howbeit",
			"however",
			"i",
			"ie",
			"if",
			"ignored",
			"immediate",
			"in",
			"inasmuch",
			"inc",
			"indeed",
			"indicate",
			"indicated",
			"indicates",
			"inner",
			"insofar",
			"instead",
			"into",
			"inward",
			"is",
			"it",
			"its",
			"itself",
			"j",
			"just",
			"k",
			"keep",
			"keeps",
			"kept",
			"know",
			"knows",
			"known",
			"l",
			"last",
			"lately",
			"later",
			"latter",
			"latterly",
			"least",
			"less",
			"lest",
			"let",
			"like",
			"liked",
			"likely",
			"little",
			"look",
			"looking",
			"looks",
			"ltd",
			"m",
			"mainly",
			"many",
			"may",
			"maybe",
			"me",
			"mean",
			"meanwhile",
			"merely",
			"might",
			"more",
			"moreover",
			"most",
			"mostly",
			"much",
			"must",
			"my",
			"myself",
			"n",
			"name",
			"namely",
			"nd",
			"near",
			"nearly",
			"necessary",
			"need",
			"needs",
			"neither",
			"never",
			"nevertheless",
			"new",
			"next",
			"nine",
			"no",
			"nobody",
			"non",
			"none",
			"noone",
			"nor",
			"normally",
			"not",
			"nothing",
			"novel",
			"now",
			"nowhere",
			"o",
			"obviously",
			"of",
			"off",
			"often",
			"oh",
			"ok",
			"okay",
			"old",
			"on",
			"once",
			"one",
			"ones",
			"only",
			"onto",
			"or",
			"other",
			"others",
			"otherwise",
			"ought",
			"our",
			"ours",
			"ourselves",
			"out",
			"outside",
			"over",
			"overall",
			"own",
			"p",
			"particular",
			"particularly",
			"per",
			"perhaps",
			"placed",
			"please",
			"plus",
			"possible",
			"presumably",
			"probably",
			"provides",
			"q",
			"que",
			"quite",
			"qv",
			"r",
			"rather",
			"rd",
			"re",
			"really",
			"reasonably",
			"regarding",
			"regardless",
			"regards",
			"relatively",
			"respectively",
			"right",
			"s",
			"said",
			"same",
			"saw",
			"say",
			"saying",
			"says",
			"second",
			"secondly",
			"see",
			"seeing",
			"seem",
			"seemed",
			"seeming",
			"seems",
			"seen",
			"self",
			"selves",
			"sensible",
			"sent",
			"serious",
			"seriously",
			"seven",
			"several",
			"shall",
			"she",
			"should",
			"since",
			"six",
			"so",
			"some",
			"somebody",
			"somehow",
			"someone",
			"something",
			"sometime",
			"sometimes",
			"somewhat",
			"somewhere",
			"soon",
			"sorry",
			"specified",
			"specify",
			"specifying",
			"still",
			"sub",
			"such",
			"sup",
			"sure",
			"t",
			"take",
			"taken",
			"tell",
			"tends",
			"th",
			"than",
			"thank",
			"thanks",
			"thanx",
			"that",
			"thats",
			"the",
			"their",
			"theirs",
			"them",
			"themselves",
			"then",
			"thence",
			"there",
			"thereafter",
			"thereby",
			"therefore",
			"therein",
			"theres",
			"thereupon",
			"these",
			"they",
			"think",
			"third",
			"this",
			"thorough",
			"thoroughly",
			"those",
			"though",
			"three",
			"through",
			"throughout",
			"thru",
			"thus",
			"to",
			"together",
			"too",
			"took",
			"toward",
			"towards",
			"tried",
			"tries",
			"truly",
			"try",
			"trying",
			"twice",
			"two",
			"u",
			"un",
			"under",
			"unfortunately",
			"unless",
			"unlikely",
			"until",
			"unto",
			"up",
			"upon",
			"us",
			"use",
			"used",
			"useful",
			"uses",
			"using",
			"usually",
			"uucp",
			"v",
			"value",
			"various",
			"very",
			"via",
			"viz",
			"vs",
			"w",
			"want",
			"wants",
			"was",
			"way",
			"we",
			"welcome",
			"well",
			"went",
			"were",
			"what",
			"whatever",
			"when",
			"whence",
			"whenever",
			"where",
			"whereafter",
			"whereas",
			"whereby",
			"wherein",
			"whereupon",
			"wherever",
			"whether",
			"which",
			"while",
			"whither",
			"who",
			"whoever",
			"whole",
			"whom",
			"whose",
			"why",
			"will",
			"willing",
			"wish",
			"with",
			"within",
			"without",
			"wonder",
			"would",
			"would",
			"x",
			"y",
			"yes",
			"yet",
			"you",
			"your",
			"yours",
			"yourself",
			"yourselves",
			"z",
			"zero",
			// stop words for paper abstracts
			//		"abstract",
			//"paper",
			//"presents",
			//"discuss",
			//"discusses",
			//"conclude",
			//"concludes",
			//"based",
			//"approach"
		};	
			//stopwords for french, added by Limin Yao
		static final String[] stopwordsFrench = {
			"fut",
			"S",
			"ces",
			"ral",
			"new",
			"tr",
			"arm",
			"y",
			"autres",
			"o",
			"tait",
			"dont",
			"ann",
			"apr",
			"sous",
			"ans",
			"cette",
			"politique",
			"of",
			"c",
			"contre",
			"leur",
			"ville",
			"fait",
			"res",
			"on",
			"deux",
			"cle",
			"v",
			"publique",
			"france",
			"te",
			"guerre",
			"sident",
			"unis",
			"mais",
			"entre",
			"aussi",
			"tat",
			"ais",
			"ses",
			"sa",
			"ont",
			"tre",
			"d",
			"pays",
			"en",
			"Il",
			"tats",
			"comme",
			"am",
			"si",
			"c",
			"fran",
			"pas",
			"g",
			"qu",
			"R",
			"aux",
			"ce",
			"f",
			"p",
			"ne",
			"son",
			"me",
			"avec",
			"l",
			"se",
			"ou",
			"sont",
			"il",
			"Les",
			"re",
			"plus",
			"m",
			"es",
			"pr",
			"la",
			"sur",
			"que",
			"pour",
			"modifier",
			"a",
			"qui",
			"Le",
			"t",
			"n",
			"au",
			"dans",
			"une",
			"par",
			"un",
			"r",
			"est",
			"e",
			"du",
			"s",
			"les",
			"en",
			"des",
			"le",
			"et",
			"l",
			"d",
			"la",
			"de",

		};

   
   
}
