<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xsl:stylesheet [
	<!ENTITY nbsp "&#160;">
	<!ENTITY deg "&#176;">
	<!ENTITY amin "&#180;">
	<!ENTITY asec "&#168;">
        <!ENTITY Beta "&#914;">
        <!ENTITY chi "&#967;">
        <!ENTITY br "&#60;&#98;&#115;&#32;&#47;&#62;">
	]>

<xsl:stylesheet version="2.0" xmlns:xsl="http://www.w3.org/1999/XSL/Transform">

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
     p.italic {font-style: italic}
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
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixLogo.png" ALT="Terapix"/>
          </TD>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixTitle.png" ALT="Logo"/>
          </TD>
          <TD ALIGN="CENTER" width="2000">
           <FONT color="#669933">
            <B> Processing summary</B>
           </FONT>
          </TD>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixPicture.gif" ALT="Terapix banner"/>
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
     <xsl:call-template name="PSFEx"/>
    </xsl:when>
   </xsl:choose>
  </xsl:for-each>
 </xsl:template>

<!-- ********************** XSL template for PSFEx *********************** -->
 <xsl:template name="PSFEx">
  <xsl:for-each select="RESOURCE[@ID='MetaData']">
   <xsl:call-template name="RunInfo"/>
   <xsl:for-each select="TABLE[@ID='PSF_Fields']">
    <xsl:call-template name="PSF_Fields"/>
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
   <b style="color: red"><xsl:if test="PARAM[@name='Error_Msg']/@value &gt; 0">
    An Error occured!!! </xsl:if>
   <xsl:value-of select="PARAM[@name='Error_Msg']/@value"/></b>
  </p>
  <p>
  <i>click to expand or hide tables</i>
  </p>
 </xsl:template>

<!-- ********************** XSL template for PSF_Fields ******************* -->
  <xsl:template name="PSF_Fields">
   <xsl:variable name="name" select="count(FIELD[@name='Catalog_Name']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ident" select="count(FIELD[@name='Image_Ident']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="next" select="count(FIELD[@name='NExtensions']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nloadtot" select="count(FIELD[@name='NStars_Loaded_Total']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="nacceptot" select="count(FIELD[@name='NStars_Accepted_Total']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fluxrad" select="count(FIELD[@name='FWHM_FromFluxRadius_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="sampling" select="count(FIELD[@name='Sampling_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="chi2" select="count(FIELD[@name='Chi2_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fwhm" select="count(FIELD[@name='FWHM_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="elongation" select="count(FIELD[@name='Elongation_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="beta" select="count(FIELD[@name='MoffatBeta_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="residuals" select="count(FIELD[@name='Residuals_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="asymmetry" select="count(FIELD[@name='Asymmetry_Mean']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="fwhmplotflag" select="count(FIELD[@name='Plot_FWHM'])"/>
   <xsl:variable name="fwhmplot" select="count(FIELD[@name='Plot_FWHM']/preceding-sibling::FIELD)+1"/>
   <xsl:variable name="ellipplotflag" select="count(FIELD[@name='Plot_Ellipticity'])"/>
   <xsl:variable name="ellipplot" select="count(FIELD[@name='Plot_Ellipticity']/preceding-sibling::FIELD)+1"/>
   <p>
    <BUTTON type="button" style="background:#CCEECC; font-family: sans-serif; font-weight: bold;" onclick="showhideTable('psfex')">
     Summary Table on Input Files
    </BUTTON>
    <TABLE class="sortable" id="psfex" BORDER="2" style="display: none">
     <TR>
      <TH BGCOLOR="#FFEECC">Filename</TH>
      <TH BGCOLOR="#FFEECC">Identifier</TH>
      <TH BGCOLOR="#FFEECC">Next</TH>
      <TH BGCOLOR="#FFEECC">Naccepted/loaded</TH>
      <TH BGCOLOR="#FFEECC">Half-Light diam.</TH>
      <TH BGCOLOR="#FFEECC">Sampling</TH>
      <TH BGCOLOR="#FFEECC">&chi;<sup>2</sup>/<i>d.o.f.</i></TH>
      <TH BGCOLOR="#FFEECC">FWHM</TH>
      <xsl:if test="$fwhmplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">FWHM map</TH>
      </xsl:if>
      <TH BGCOLOR="#FFEECC">Elongation</TH>
      <xsl:if test="$ellipplotflag &gt; 0">
       <TH BGCOLOR="#FFEECC">Ellipticity map</TH>
      </xsl:if> 
      <TH BGCOLOR="#FFEECC">&Beta;</TH>
      <TH BGCOLOR="#FFEECC">Residuals</TH>
      <TH BGCOLOR="#FFEECC">Asymmetry</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
<!-- Catalog name -->
        <td  BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$name]"/></el>
        </td>
<!-- Identifier -->
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$ident]"/></el>
        </td>
<!-- Number of extensions -->
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$next]"/></el>
        </td>
<!-- Number of sources -->
        <td align="center" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="TD[$nacceptot]"/></el> / <el><xsl:value-of select="TD[$nloadtot]"/></el>
        </td>
<!-- Half-light diameter -->
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$fluxrad],'##0.000')"/></el>
        </td>
<!-- Sampling -->
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$sampling],'##0.000')"/></el>
        </td>
<!-- Chi2 --> 
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$chi2], '######0.0')"/></el>
        </td>
<!-- FWHM -->
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$fwhm],'##0.000')"/></el>
        </td>
<!-- FWHM map -->
        <xsl:if test="$fwhmplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$fwhmplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$fwhmplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$fwhmplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
<!-- Elongation -->
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$elongation],'##0.000')"/></el>
        </td>
<!-- Ellipticity map -->
        <xsl:if test="$ellipplotflag &gt; 0">
         <td align="center" BGCOLOR="#EEEEEE">
          <a target="_blank">
           <xsl:attribute name="href">
            <xsl:value-of select="TD[$ellipplot]"/>
           </xsl:attribute>
           <img width="80">
            <xsl:attribute name="src">
             <xsl:value-of select="TD[$ellipplot]"/>
            </xsl:attribute>
            <xsl:attribute name="title">
             <xsl:value-of select="TD[$ellipplot]"/>
            </xsl:attribute>
           </img>
          </a>
         </td>
        </xsl:if>
<!-- Beta -->
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$beta],'##0.000')"/></el>
        </td>
<!-- Residuals -->
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$residuals],'##0.000')"/></el>
        </td>
<!-- Asymmetry -->
        <td align="right" BGCOLOR="#EEEEEE">
         <el><xsl:value-of select="format-number(TD[$asymmetry],'##0.000')"/></el>
        </td>
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
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
     Warnings (limited to the last 1000)
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

