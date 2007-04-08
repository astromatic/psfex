<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xsl:stylesheet [
	<!ENTITY nbsp "&#160;">
	<!ENTITY deg "&#176;">
	<!ENTITY amin "&#180;">
	<!ENTITY asec "&#168;">
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
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixLogo.png" ALT="Terapix"/>
          </TD>
          <TD ALIGN="CENTER">
           <IMG SRC="http://terapix.iap.fr/cplt/xsl/terapixTitle.png" ALT="Logo"/>
          </TD>
          <TD ALIGN="CENTER">
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
 </xsl:template>

<!-- ********************** XSL template for PSF ********************** -->
  <xsl:template name="PSF">
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
   <p>
    <TABLE class="sortable" id="psfex" BORDER="2">
     <TR>
      <TH BGCOLOR="#FFEECC">N. star loaded</TH>
      <TH BGCOLOR="#FFEECC">N. star accepted</TH>
      <TH BGCOLOR="#FFEECC">FWHM from Flux Radius</TH>
      <TH BGCOLOR="#FFEECC">Sampling</TH>
      <TH BGCOLOR="#FFEECC">Chi2</TH>
      <TH BGCOLOR="#FFEECC">Min FWHM</TH>
      <TH BGCOLOR="#FFEECC">FWHM</TH>
      <TH BGCOLOR="#FFEECC">Max FWHM</TH>
      <TH BGCOLOR="#FFEECC">Min Elongation</TH>
      <TH BGCOLOR="#FFEECC">Elongation</TH>
      <TH BGCOLOR="#FFEECC">Max Elongation</TH>
      <TH BGCOLOR="#FFEECC">Min Beta</TH>
      <TH BGCOLOR="#FFEECC">Beta exponent</TH>
      <TH BGCOLOR="#FFEECC">Max Beta</TH>
      <TH BGCOLOR="#FFEECC">Min Residuals</TH>
      <TH BGCOLOR="#FFEECC">Residuals</TH>
      <TH BGCOLOR="#FFEECC">Max Residuals</TH>
     </TR>
     <xsl:for-each select="DATA/TABLEDATA">
      <xsl:for-each select="TR">
       <tr>
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
       </tr>
      </xsl:for-each>
     </xsl:for-each>
    </TABLE>
   </p>
 </xsl:template>


 <xsl:template name="Rest">
</xsl:template>
</xsl:stylesheet>
