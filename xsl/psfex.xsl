<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xsl:stylesheet [
	<!ENTITY nbsp "&#160;">
	<!ENTITY deg "&#176;">
	<!ENTITY amin "&#180;">
	<!ENTITY asec "&#168;">
        <!ENTITY chi "&#967;">
        <!ENTITY beta "&#946;">
        <!ENTITY copy "&#169;">
        <!ENTITY br "&#xD;&#xA;">
	]>

<xsl:stylesheet version="1.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

<!-- *********************** Global XSL template ************************** -->
 <xsl:template match="/">
  <xsl:variable name="date" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Date']/@value"/>
  <xsl:variable name="time" select="/VOTABLE/RESOURCE/RESOURCE[@name='MetaData']/PARAM[@name='Time']/@value"/>
  <html>

<!-- HTML head  -->

   <head>

<!-- javascript -->

<!--  <script type="text/javascript" language="javascript"> -->

    <script src="http://terapix.iap.fr/cplt/xsl/sorttable.js"/>

    <style type="text/css">
     p.sansserif {font-family: sans-serif}
     body {background-color: white}
     mono {font-family: monospace}
     elen {font-family: monospace; font-size: 100%; font-weight: bold; color: green }
     elep {font-family: monospace; font-size: 100%; font-weight: bold; color: red }
     el {font-family: monospace; font-size: 100%; color: black}
     a {text-decoration: none}
     table.sortable a.sortheader
      {
      background-color:#FFEECC;
      color: black;
      font-weight: bold;
      font-size: 80%;
      text-decoration: none;
      display: button;
      }
     table.sortable span.sortarrow
      {
      color: black;
      font-weight: bold;
      text-decoration: none;
      }
     table.sortable a.sortheader.sub
      {
      vertical-align: sub;
      }
     </style>

     <title>
      Processing summary on <xsl:value-of select="$date"/> at <xsl:value-of select="$time"/>
     </title>
    </head>

<!-- HTML body -->

    <BODY>
     <TABLE BORDER="0" CELLPADDING="0" CELLSPACING="0" WIDTH="100%">
      <TR>
       <TD ALIGN="LEFT">
        <TABLE BORDER="0">
         <TR>
          <TD ALIGN="CENTER">
           <a href="http://terapix.iap.fr"><IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixLogo.png" ALT="Terapix" title="Terapix" border="0"/></a>
          </TD>
          <TD ALIGN="CENTER">
           <a href="http://terapix.iap.fr"><IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixTitle.png" ALT="Logo" title="Terapix" border="0"/></a>
          </TD>
          <TD ALIGN="CENTER" width="2000">
           <FONT color="#669933">
            <B> Processing summary</B>
           </FONT>
          </TD>
          <TD ALIGN="CENTER">
           <IMG
     SRC="http://terapix.iap.fr/cplt/xsl/terapixPicture.gif" title="Star Formation Region IC 1396, &copy; 2002 CFHT" alt="Star Formation Region IC 1396, &copy; 2002 CFHT"/>
          </TD>
         </TR>
        </TABLE>
       </TD>
      </TR>
      <TR>
       <TD>
        <TABLE BORDER="0" WIDTH="100%" BGCOLOR="#000000">
         <TR>
          <TH BGCOLOR="#000000" ALIGN="LEFT"><FONT SIZE="-1" COLOR="#FFFFFF"> Home > Tools > Data reduction</FONT></TH>
         </TR>
        </TABLE>
       </TD>
      </TR>
     </TABLE>
    <xsl:call-template name="VOTable"/>
   </BODY>
  </html>
 </xsl:template>

<!-- **************** Generic XSL template for VOTables ****************** -->
 <xsl:template name="VOTable">
  <xsl:for-each select="/VOTABLE">
   <xsl:call-template name="Resource"/>
  </xsl:for-each>
 </xsl:template>

<!-- *************** Generic XSL template for Resources ****************** -->
 <xsl:template name="Resource">
  <xsl:for-each select="RESOURCE">
   <xsl:choose>
    <xsl:when test="@ID='PSFEx'">
     <xsl:call-template name="psfex"/>
    </xsl:when>
   </xsl:choose>
  </xsl:for-each>
 </xsl:template>

<!-- ********************** XSL template for PSFEx *********************** -->
 <xsl:template name="psfex">
  <xsl:for-each select="RESOURCE[@ID='MetaData']">
   <xsl:call-template name="RunInfo"/>
   <xsl:for-each select="TABLE[@ID='PSF']">
    <xsl:call-template name="PSF"/>
   </xsl:for-each>
   <xsl:for-each select="RESOURCE[@ID='Config']">
    <xsl:call-template name="Config"/>
   </xsl:for-each>
   <xsl:for-each select="TABLE[@ID='Warnings']">
    <xsl:call-template name="Warnings"/>
   </xsl:for-each>
  </xsl:for-each>
 </xsl:template>

<!-- ************* Generic XSL RunInfo template for MetaData ************* -->
 <xsl:template name="RunInfo">
  <p>
<!-- Software name, version, date, time and number of threads -->
   <a>
    <xsl:attribute name="href">
     <xsl:value-of select="PARAM[@name='Soft_URL']/@value"/>
    </xsl:attribute>
    <b>
     <xsl:value-of select="PARAM[@name='Software']/@value"/>&nbsp;<xsl:value-of select="PARAM[@name='Version']/@value"/>
    </b>
   </a>
   started on
   <b><xsl:value-of select="PARAM[@name='Date']/@value"/></b>
   at
   <b><xsl:value-of select="PARAM[@name='Time']/@value"/></b>
   with
   <b><xsl:value-of select="PARAM[@name='NThreads']/@value"/></b>
   thread<xsl:if test="PARAM[@name='NThreads']/@value &gt; 1">s</xsl:if>

<!-- Run time -->
   <xsl:variable name="duration" select="PARAM[@name='Duration']/@value"/>
   (run time:
    <b>
     <xsl:choose> 
      <xsl:when test="$duration &gt; 3600.0">
       <xsl:value-of
	select='concat(string(floor($duration div 3600)),
	" h ", format-number(floor(($duration div 60) mod 60.0), "00"),
	" min")'/>
      </xsl:when>
      <xsl:otherwise>
       <xsl:choose>
        <xsl:when test="$duration &gt; 60.0">
         <xsl:value-of
	  select='concat(format-number(floor($duration div 60),"##"),
	  " min ", format-number(floor($duration mod 60.0), "00")," s")'/>
        </xsl:when>
        <xsl:otherwise>
         <xsl:value-of select='concat(string($duration), " s")'/>
        </xsl:otherwise>
       </xsl:choose>
      </xsl:otherwise>
     </xsl:choose>
    </b>)
    <br />
   by user <b><xsl:value-of select="PARAM[@name='User']/@value"/></b>
   from <b><xsl:value-of select="PARAM[@name='Host']/@value"/></b>
   in <b><mono><xsl:value-of select="PARAM[@name='Path']/@value"/></mono></b>
  </p>
  <p>
  <sans-serif><i>click to expand or hide tables</i></sans-serif>
  </p>
 </xsl:template>

<!-- ********************** XSL template for PSF ********************** -->
  <xsl:template name="PSF">
   <xsl:variable name="ext" select="count(FIELD[@name='Extension']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nload" select="count(FIELD[@name='NStars_Loaded']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="naccept" select="count(FIELD[@name='NStars_Accepted']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fwhmfr" select="count(FIELD[@name='FWHM_FromFluxRadius']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="sampling" select="count(FIELD[@name='Sampling']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="chi2" select="count(FIELD[@name='Chi2']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fwhm_min" select="count(FIELD[@name='FWHM_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fwhm" select="count(FIELD[@name='FWHM']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fwhm_max" select="count(FIELD[@name='FWHM_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="elong_min" select="count(FIELD[@name='Elongation_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="elong" select="count(FIELD[@name='Elongation']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="elong_max" select="count(FIELD[@name='Elongation_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="beta_min" select="count(FIELD[@name='MoffatBeta_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="beta" select="count(FIELD[@name='MoffatBeta']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="beta_max" select="count(FIELD[@name='MoffatBeta_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="res_min" select="count(FIELD[@name='Residuals_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="res" select="count(FIELD[@name='Residuals']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="res_max" select="count(FIELD[@name='Residuals_Max']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="as_min" select="count(FIELD[@name='Asymmetry_Min']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="as" select="count(FIELD[@name='Asymmetry']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="as_max" select="count(FIELD[@name='Asymmetry_Max']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: sans-serif; font-weight: bold;" onclick="showhideTable('psfex')">
     Summary Table on PSF Model
    </BUTTON>
    <TABLE class="sortable" id="psfex" BORDER="2" style="display: none">
     <TR>
      <TH BGCOLOR="#FFEECC">Extension</TH>
      <TH BGCOLOR="#FFEECC">N. star loaded</TH>
      <TH BGCOLOR="#FFEECC">N. star accepted</TH>
      <TH BGCOLOR="#FFEECC">FWHM from Flux Radius &br; (<xsl:value-of select="FIELD[@name='FWHM_FromFluxRadius']/@unit"/>)</TH>
      <TH BGCOLOR="#FFEECC">Sampling &br; (<xsl:value-of select="FIELD[@name='Sampling']/@unit"/>)</TH>
      <TH BGCOLOR="#FFEECC">&chi;<sup>2</sup></TH>
      <TH BGCOLOR="#FFEECC">FWHM min &br; (<xsl:value-of select="FIELD[@name='FWHM_Min']/@unit"/>)</TH>
      <TH BGCOLOR="#FFEECC">FWHM mean &br; (<xsl:value-of select="FIELD[@name='FWHM']/@unit"/>)</TH>
      <TH BGCOLOR="#FFEECC">FWHM max &br; (<xsl:value-of select="FIELD[@name='FWHM_Max']/@unit"/>)</TH>
      <TH BGCOLOR="#FFEECC">Elongation min</TH>
      <TH BGCOLOR="#FFEECC">Elongation mean</TH>
      <TH BGCOLOR="#FFEECC">Elongation max</TH>
      <TH BGCOLOR="#FFEECC">Moffat &beta; min</TH>
      <TH BGCOLOR="#FFEECC">Moffat &beta; mean</TH>
      <TH BGCOLOR="#FFEECC">Moffat &beta; max</TH>
      <TH BGCOLOR="#FFEECC">Residuals min</TH>
      <TH BGCOLOR="#FFEECC">Residuals mean</TH>
      <TH BGCOLOR="#FFEECC">Residuals max</TH>
      <TH BGCOLOR="#FFEECC">Asymmetry min</TH>
      <TH BGCOLOR="#FFEECC">Asymmetry mean</TH>
      <TH BGCOLOR="#FFEECC">Asymmetry max</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$ext]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$nload]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$naccept]"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$fwhmfr], '#0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$sampling], '#0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$chi2], '#0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$fwhm_min], '#0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$fwhm], '#0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$fwhm_max], '#0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$elong_min], '#0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$elong], '#0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$elong_max], '#0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$beta_min], '#0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$beta], '#0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$beta_max], '#0.000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$res_min], '#0.0000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$res], '#0.0000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$res_max], '#0.0000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$as_min], '#0.0000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$as], '#0.0000')"/></el>
        </td>
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$as_max], '#0.0000')"/></el>
        </td>
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
   <xsl:if test="PARAM[@name='NExtensions']/@value &gt; 1">
    <p>
     <BUTTON type="button" style="background:#CCEECC; font-family: sans-serif; font-weight: bold;" onclick="showhideTable('avpsfex')">
      Average Properties of PSF on all Extensions
     </BUTTON>
     <TABLE id="avpsfex" class="sortable" BORDER="2" style="display: none">
      <TR>
       <TH BGCOLOR="#FFEECC"> </TH>
       <TH BGCOLOR="#FFEECC">N. star loaded</TH>
       <TH BGCOLOR="#FFEECC">N. star accepted</TH>
       <TH BGCOLOR="#FFEECC">FWHM from Flux Radius &br; (<xsl:value-of select="PARAM[@name='FWHM_FromFluxRadius_Mean']/@unit"/>)</TH>
       <TH BGCOLOR="#FFEECC">Sampling &br; (<xsl:value-of select="PARAM[@name='Sampling_Mean']/@unit"/>)</TH>
       <TH BGCOLOR="#FFEECC">&chi;<sup>2</sup></TH>
       <TH BGCOLOR="#FFEECC">FWHM &br; (<xsl:value-of select="PARAM[@name='FWHM_Mean']/@unit"/>)</TH>
       <TH BGCOLOR="#FFEECC">Elongation</TH>
       <TH BGCOLOR="#FFEECC">Moffat &beta;</TH>
       <TH BGCOLOR="#FFEECC">Residuals</TH>
       <TH BGCOLOR="#FFEECC">Asymmetry</TH>
      </TR>
      <tr BGCOLOR="#EEEEEE">
       <td BGCOLOR="#FFEECC"><b>min</b></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='NStars_Loaded_Min']/@value,'#0.00')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='NStars_Accepted_Min']/@value,'#0.00')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='FWHM_FromFluxRadius_Min']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Sampling_Min']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Chi2_Min']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='FWHM_Min']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Elongation_Min']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='MoffatBeta_Min']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Residuals_Min']/@value,'#0.0000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Asymmetry_Min']/@value,'#0.0000')"/></el></td>
      </tr>
      <tr BGCOLOR="#EEEEEE">
       <td BGCOLOR="#FFEECC"><b>mean</b></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='NStars_Loaded_Mean']/@value,'#0.00')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='NStars_Accepted_Mean']/@value,'#0.00')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='FWHM_FromFluxRadius_Mean']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Sampling_Mean']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Chi2_Mean']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='FWHM_Mean']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Elongation_Mean']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='MoffatBeta_Mean']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Residuals_Mean']/@value,'#0.0000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Asymmetry_Mean']/@value,'#0.0000')"/></el></td>
      </tr>
      <tr BGCOLOR="#EEEEEE">
       <td BGCOLOR="#FFEECC"><b>max</b></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='NStars_Loaded_Max']/@value,'#0.00')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='NStars_Accepted_Max']/@value,'#0.00')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='FWHM_FromFluxRadius_Max']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Sampling_Max']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Chi2_Max']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='FWHM_Max']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Elongation_Max']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='MoffatBeta_Max']/@value,'#0.000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Residuals_Max']/@value,'#0.0000')"/></el></td>
       <td align="right"><el><xsl:value-of select="format-number(PARAM[@name='Asymmetry_Max']/@value,'#0.0000')"/></el></td>
      </tr>
     </TABLE>
    </p>
   </xsl:if>
 </xsl:template>

<!-- ********************** XSL template for Config File ********************** -->
  <xsl:template name="Config">
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: sans-serif; font-weight: bold;" onclick="showhideTable('config')">
     Configuration File: <xsl:value-of select="PARAM[@name='Prefs_Name']/@value"/>
    </BUTTON>
    <TABLE id="config" class="sortable" style="display: none">
     <TR>
      <TH BGCOLOR="#FFEECC">Config Parameter</TH>
      <TH BGCOLOR="#FFEECC">Value</TH>
     </TR>
     <xsl:for-each select="PARAM[position()>2]">
      <tr BGCOLOR="#EEEEEE">
       <td><el><xsl:value-of select="@name"/></el></td>
       <td><el><xsl:value-of select="@value"/></el></td>
      </tr>
     </xsl:for-each>
    </TABLE>
   </p>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: monospace; font-weight: bold: font-size: 80%;" onclick="showhideTable('commandline')">
     Command Line
    </BUTTON>
    <TABLE id="commandline" style="display: none">
     <TR>
      <TD BGCOLOR="#FFEECC" style="font-size: 80%;"><el>Command Line: <xsl:value-of select="PARAM[@name='Command_Line']/@value"/></el></TD>
     </TR>
    </TABLE>
   </p>
  </xsl:template>

<!-- ********************** XSL template for Warnings ********************** -->
  <xsl:template name="Warnings">
   <xsl:variable name="date" select="count(FIELD[@name='Date']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="time" select="count(FIELD[@name='Time']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="msg" select="count(FIELD[@name='Msg']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: monospace; font-weight: bold: font-size: 80%;" onclick="showhideTable('warnings')">
     Warnings (limited to the last 100)
    </BUTTON>
    <TABLE id="warnings" style="display: none">
     <TR style="font-size: 80%;">
      <TH BGCOLOR="#FFEECC">Date</TH>
      <TH BGCOLOR="#FFEECC">Time</TH>
      <TH BGCOLOR="#FFEECC">Message</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$date]"/></el>
        </td>
        <td BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$time]"/></el>
        </td>
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$msg]"/></el>
        </td>
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
 </xsl:template>


 <xsl:template name="Rest">
</xsl:template>
</xsl:stylesheet>
