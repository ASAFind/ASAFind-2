
### ![ASAFind logo](ASAFind_logo_small.png "ASAFind logo")

<font color="#00c03a"> <h3>Download ASAFind to install locally</h3></font>

<p>For local installation, a command line version of ASAFind can be downloaded from our GitHub repository:</p>
<a href="https://github.com/ASAFind/ASAFind-2" target="_blank"> https://github.com/ASAFind/ASAFind-2</a>
<p>
<font color="#00c03a"> <h5>Installation steps</h5></font></p>
Using the GitHub function, you can either dowload the files as a zip archive, or you can clone the repository using the provided URL. After download of the latest version, follow these installation steps:

<ul>
    <li>step into the directory where you want to make the installation e.g. <br> <font color="#00c03a"> cd /home/marta/asafind</font></li>
    <li>Make a clone of the GitHub  <br><font color="#00c03a"> git clone https://github.com/ASAFind/ASAFind-2.git </font></li>
    <li>run the following command from the command line <br><font color="#00c03a"> python3 -m venv asafind_line_command 
     . asafind_line_command/bin/activate <br>
     pip install --upgrade pip <br>
     cd ASAFind-2/ <br>
     pip install -r requirements.txt <br></font></li>
    <li>Now you are in virtual environment named asafind_line_command. Here you can ask for help e.g.
    <br><font color="#00c03a">(asafind_line_command) marta@marta-mini:~/asafind/ASAFind-2$ python3 S0_ASAFind.py  --help</font></li>

</ul>

It will create the environment called asafind_line_command including subdirectries temp and output, install all required packages activate it.

### ![ASAFind structure](directories.png "ASAFind structure")

<br>

<font color="#00c03a"> <h4>ASAFind 2.0</h4></font><br>
In the environment root directory is the program S1_ASAFind_v3.py<br>
Takes a Fasta and companion TargetP v.2.0 short format file as input, with the complete TargetP header <br>
(two lines starting with '#'). Some versions of SignalP truncate the sequence names. SignalP-3.0 to 20 characters, <br>
and 4.0, 4.1 to 58 characters. Therefore,<br>
ASAFind only considers the first corresponding characters of the fasta name (and the first 90 in the<br>
case of TargetP 2.0), which must be unique within the file. Parts of the fasta name<br>
after that character are ignored. Additionally, the fasta name may not contain a '-' or '|'. This<br>
requirement is because SignalP converts characters in sequence names (e.g. '-' is changed to '_').<br>
ASAFind requires at least 7 aa upstream and 22 aa downstream of the cleavage site suggested by<br>
SignalP. The output of this script is a tab delimited table.<br>
Python >= 3.10 required.<br><br>

<font color="#52595D"> 
<b>python S0_ASAFind.py --help</b>
<br>

<table style="width:100%; border-collapse: collapse;">
    <tr>
      <td style="border: 1px solid black;">usage: S1_ASAFind_v3.py</td>
      <td style="border: 1px solid black;"></td>
    </tr>
    <tr>
      <td style="border: 1px solid black;">usage: S0_ASAFind.py</td>
      <td style="border: 1px solid black;">
        [-h] -f FASTA_FILE -p SIGNALP_FILE<br>
        [-s SIMPLE_SCORE_CUTOFF] [-t FASTA_FILE_WITH_MOTIFS] [-w]<br>
        [-v1] [-ppc] [-s_ppc SCORE_CUTOFF_PPC]<br>
        [-t_ppc FASTA_FILE_WITH_MOTIFS_PPC] [-l]<br>
        [-my_org MY_ORGANISM] [-v]
      </td>
    </tr>
  </table>
  
<br><br>

<table style="width:100%; border-collapse: collapse;">
<tr>
  <td style="border: 1px solid black;"> -h, --help</td><td style="border: 1px solid black;">show this help message and exit</td>
</tr><tr>
  <td style="border: 1px solid black;">-f FASTA_FILE, --fasta_file FASTA_FILE</td><td style="border: 1px solid black;">
                        Specify the input fasta FILE.</td>
</tr><tr>                    
   <td style="border: 1px solid black;">-p SIGNALP_FILE, --signalp_file SIGNALP_FILE</td><td style="border: 1px solid black;">
                        Specify the input TargetP FILE..</td>
</tr><tr>
   <td style="border: 1px solid black;">-s SIMPLE_SCORE_CUTOFF, --simple_score_cutoff SIMPLE_SCORE_CUTOFF</td><td style="border: 1px solid black;">
                        Optionally, specify an explicit score cutoff, rather
                        than using ASAFind's default algorithm, not compatible
                        with option -v1. The score given here will not be
                        normalized and therefore should be obtained form a
                        distribution of normalized scores.</td>
</tr><tr>
    <td style="border: 1px solid black;">-t FASTA_FILE_WITH_MOTIFS, --fasta_file_with_motifs FASTA_FILE_WITH_MOTIFS</td><td style="border: 1px solid black;">
                        Optionally, specify a custom scoring table. The
                        scoring table will be normalized with the maximum
                        score, which allows for processing of non-normalized
                        as well as normalized scoring tables.</td>
</tr><tr>
    <td style="border: 1px solid black;">-w, --web_output</td><td style="border: 1px solid black;">
                        Format output for web display. This is mostly useful
                        when called by a web app. </td>                       

</tr><tr>
    <td style="border: 1px solid black;">-v1, --reproduce_ASAFind_1</td><td style="border: 1px solid black;">
                        Reproduce ASAFind 1.x scores and results (non-normalized scores, if no custom scoring table is<br>
                        specified, the original default scoring table generated without small sample size correction<br>
                        will be used, not compatible with option -s).</td>

</tr><tr>
    <td style="border: 1px solid black;">-ppc, --include_ppc_prediction</td><td style="border: 1px solid black;">
                        Include prediction of proteins that might be targeted to the periplastidic compartment.</td>


</tr><tr>                    
    <td style="border: 1px solid black;">-t SCORE_TABLE_FILE, --score_table_file SCORE_TABLE_FILE</td>
    <td style="border: 1px solid black;">Optionally, specify a custom scoring table. The scoring table will be normalized with the<br>
                        maximum score, which allows for processing of non-normalized as well as normalized scoring<br>
                        tables.</td>
</tr><tr>
    <td style="border: 1px solid black;">-o OUT_FILE, --out_file OUT_FILE</td>
    <td style="border: 1px solid black;">Specify the path and name of the output file you wish to create. Default will be the same as<br>
                        the fasta_file, but with a ".tab" suffix.</td>

</tr><tr>
    <td style="border: 1px solid black;">-s_ppc SCORE_CUTOFF_PPC, --score_cutoff_ppc SCORE_CUTOFF_PPC</td>
    <td style="border: 1px solid black;">Optionally, specify an explicit score cutoff for the ppc protein prediction, if given, ppc<br>
                        protein prediction will be included. The score given here will not be normalized and therefore<br>
                        should be obtained form a distribution of normalized scores.</td>

</tr><tr>
    <td style="border: 1px solid black;">-t_ppc SCORE_TABLE_FILE_PPC, --score_table_file_ppc SCORE_TABLE_FILE_PPC</td>
    <td style="border: 1px solid black;">Optionally, specify a custom scoring table for the ppc protein prediction, if given, ppc<br>
                        protein prediction will be included. The scoring table will be normalized with the maximum<br>
                        score, which allows for processing of non-normalized as well as normalized scoring tables.</td>
</tr><tr>
    <td style="border: 1px solid black;">-l {yes,no}, --logomaker {yes,no}</td>
    <td style="border: 1px solid black;">Choose "yes" or "no" (default: no). Optionally, from version 3.0 you can define with keyword<br>
                        yes or no, if the program will generate also the logomaker pictures in .png and .svg formats.<br>
                        They will be include into the output compressed package.</td>
</tr><tr>
    <td style="border: 1px solid black;">-my_org MY_ORGANISM, --my_organism MY_ORGANISM</td>
    <td style="border: 1px solid black;">Specify the name of organism.</td>
</tr><tr>
    <td style="border: 1px solid black;">-v, --version </td>        <td style="border: 1px solid black;">show program's version number and exit<td>
    </tr>
</table>
<br>
<font color="#00c03a"> <h4>Example of run</h4>

python S0_ASAFind.py -f /data/haptophyta/temp/haptophyta_fasta_for_targetp2.fsa -p /data/haptophyta/temp/haptophyta_fasta_for_targetp2.targetp2 -my_org haptophyta -l <br><br>
</font>
