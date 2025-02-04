#/data/software/anaconda3/envs/lhri_bioinformatics/bin/ python

import sys
import os
from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Frame, Image, PageBreak, PageTemplate, Paragraph, Spacer, Table, TableStyle, XPreformatted, KeepTogether
from sierrapy import fastareader, SierraClient
from sierrapy.sierraclient import ResponseError
from textwrap import TextWrapper
from datetime import datetime
from gql import gql
import json
from helper import *
from Crimson_metadata import *
import re
from Bio import SeqIO

class ConditionalSpacer(Spacer):
    def wrap(self, availWidth, availHeight):
        height = min(self.height, availHeight-1e-8)
        return (availWidth, height)

def main(fasta, response):
    #Makes API call to stanford's hivdb
    sequences = {}
    dir_path = os.path.dirname(os.path.realpath(__file__))
    for seq_record in SeqIO.parse(fasta, "fasta"):
        sequences['header']=seq_record.id
        sequences['sequence']=(str(seq_record.seq))
    #print(sequences['sequence'])
    #sequences['header'], sequences['sequence'] =open(fasta).decode("ascii")[1:].splitlines()
    #sequences['header'], sequences['sequence'] =fastareader.load(open(fasta))[1:].splitlines()
    client = SierraClient('http://localhost:8080/sierra/rest/graphql')
    # client = SierraClient('https://hivdb.stanford.edu/graphql')
    # try:
    fileIn = client.execute(gql(open(os.path.join(dir_path,'default.gql')).read()),variable_values={"sequences": [sequences]})

    seqId=sequences['header'].split(".")[0]
    if len(re.findall(r'P\d{6}', seqId)) > 0:
        sampleId = re.findall(r'P\d{6}', seqId)[0]
    else:
        sampleId = seqId

    metadata = getMetadata(sampleId)

    expectedKeys = ['Patient Name','NIH ID', 'Sample ID','pnum_noP','Sample Date']
    for i,val in enumerate(expectedKeys):
        if not val in metadata:
            intKey="_"+str(i)
            if intKey in metadata:
                metadata[val] = metadata[intKey]
                metadata.pop(intKey, None)

    metadata["Sample ID"] = sampleId
    metadata["Sample Type"] = "Plasma"
    metadata["Sequence ID"] = seqId
    #Mark sequences variable to be garbage collected
    del sequences
    version = fileIn["viewer"]["currentVersion"]['text']
    json_out=os.path.splitext(response)[0]+'.json'

    json.dump(fileIn, open(json_out,'w'),separators=(',', ': '), indent = 2)

# #########################################################################################
# Added this block on 2/20/2020 to add the patient name to the pdf file name
    patient_name = metadata.get('Patient Name')
    if patient_name and patient_name.strip():
        patient_name_list=patient_name.split(",")
        response_list=[os.path.dirname(response),os.path.basename(response)]
        response=response_list[0]+"/"+patient_name_list[0]+"_"+patient_name_list[1][0:1]+"-"+response_list[1]
        # print("\n*******************"+response_list[0]+"/"+patient_name_list[0]+"_"+patient_name_list[1][0:1]+"_"+response_list[1]+"**********************\n")
    # else:
    #     print("\n*******************"+response+"**********************\n")

# #########################################################################################


    #Remove extra data from dict that is no longer needed
    fileIn = fileIn["viewer"]["sequenceAnalysis"][0]
    headerData = {}
    #fileIn = json.load(open("./reports/Positive-control-RUO_S3.Amb5.json","r"))
    #==============================================================================================================================================
    #REFORMAT JSON DATA TO CREATE TABLES
    #==============================================================================================================================================

    nn = True #placeholder toggle
    #dictionary of 2d arrays to store data for reportlab tables
    tableData = {}
    for dict in fileIn["drugResistance"]: #Each dictionary in drugResistance is a separate drug Class
        gene = dict["gene"]["name"]
        if gene == "RT":#Splits "RT" dict into two subdicts for NNRTI and NRTI
            tableData["NRTI"] = []
            tableData["NRTI"].append(["Drug", "Mutations List", "Score", "Range","Color", "Interpretation"])
        else: #Standard behavior for IN and PR
            tableData[gene] = []
            tableData[gene].append(["Drug", "Mutations List", "Score", "Range", "Color", "Interpretation"])
        for drug in dict["drugScores"]:#Iterates Through each drug in class
            mutations = set([]) #Set of unique mutations
            for m in drug["partialScores"]:
                for mut in m["mutations"]:
                    mutations.add(mut["text"]) #Gets set of UNIQUE mutations that affect susceptibility
            mutations = list(mutations)
            #print(mutations)
            #Sorts mutations by location
            for i in range(len(mutations)-1):
                for j in range(len(mutations)-1):
                    if int(list(filter(str.isdigit, mutations[j]))[0]) > int(list(filter(str.isdigit, mutations[j+1]))[0]):
                        mutations[j], mutations[j+1] = mutations[j+1], mutations[j]
            mutations = ", ".join(mutations)
            #Separates NRTI and NNRTI
            if nn and drug["drugClass"]["name"] == "NNRTI":
                nn=False
                tableData["NNRTI"] = []
                tableData["NNRTI"].append(["Drug", "Mutations List", "Score", "Range", "Color", "Interpretation"])
            if gene == "RT":
                if not nn:
                    tableData["NNRTI"].append([drug["drug"]["fullName"]+" ("+drug["drug"]["displayAbbr"]+")",
                        mutations, int(drug["score"]), int(drug["level"]), "", drug["text"]])
                else:
                    tableData["NRTI"].append([drug["drug"]["fullName"]+" ("+drug["drug"]["displayAbbr"]+")",
                        mutations, int(drug["score"]), int(drug["level"]), "", drug["text"]])
            else:
                tableData[gene].append([drug["drug"]["fullName"]+" ("+drug["drug"]["displayAbbr"]+")",
                    mutations, int(drug["score"]), int(drug["level"]), "", drug["text"]])
    #==============================================================================================================================================
    #CREATE PDF
    #==============================================================================================================================================
    #Register times new roman font
    pdfmetrics.registerFont(TTFont('Times',os.path.join(dir_path,'times.ttf')))
    pdfmetrics.registerFont(TTFont('TimesBd', os.path.join(dir_path,'timesbd.ttf')))
    pdfmetrics.registerFontFamily('Times', normal='Times', bold='TimesBd',
                                  italic='Times', boldItalic='TimesBd')

    #Additional stuff shown in footer of first page
    def _header_footerFP(canvas, doc):
        pos = _header_footer(canvas, doc)#Gets the coordinate of the top of the footer
        canvas.saveState()
        #Fixed categories for the table
        fpTableCats = [["Signatures:",""],
                       ["Performed By:","Reviewed by:"],
                       ["Report Date:", "Review Date:"]]
        fpTableVal = [["",""],["MK/GM",""],[datetime.now().strftime("%m/%d/%y"),""]]
        #Converts categories and data to formatted paragraphs
        fpTableData = []
        for i, vals in enumerate(fpTableCats):
            fpTableData.append([])
            for j, val in enumerate(vals):
                fpTableData[i].append(Paragraph('<b>'+val+' </b>'+fpTableVal[i][j], style=styles['std9']))
        #Style for the table
        tblStyle = TableStyle([('FONTSIZE',(0,0),(-1,-1), 8),
                               ('FONT', (0,0), (-1,-1), 'Times'),
                               ('VALIGN',(0,0),(-1,-1),'TOP'),
                               ('BOTTOMPADDING',(0,0),(-1,-1),0),
                               ('TOPPADDING',(0,0),(-1,-1),0),
                               ('LEFTPADDING',(0,0),(-1,-1),0),
                               ('ENDPADDING',(0,0),(-1,-1),0)])
        #Creates the table
        fpTable = Table(fpTableData,colWidths=[396,144],rowHeights=[10.8,10.8,10.8],
            style=tblStyle)
        fpTable.wrap(doc.width, doc.topMargin)
        #Draws table on canvas
        fpTable.drawOn(canvas, doc.leftMargin,pos)

        #Loads and draws signature
        sigs = Image(os.path.join(dir_path,"RobinSignature.jpg"), width=90,height=22.7)
        sigs.wrap(doc.width, doc.topMargin)
        sigs.drawOn(canvas, 90, 85)
        canvas.restoreState()


    def _header_footer(canvas, doc):
        canvas.saveState()
        #HEADER
        #Top Banner Image
        hImg = Image(os.path.join(dir_path,"header.bmp"), width=doc.width, height = 59.7)
        imgW, imgH = hImg.wrap(doc.width, doc.topMargin)
        hImg.drawOn(canvas, doc.leftMargin, doc.height - doc.topMargin - imgH)

        #First Horizontal Rule
        hL1 = lineFlow(doc._rightMargin - doc.leftMargin, 0)
        lW1, lH1 = hL1.wrap(doc.width, doc.topMargin) #Gets Width and height of line
        #Gets Vertical coordinate of line
        L1Pos = doc.height - doc.topMargin - lH1 - imgH - 10.4
        hL1.drawOn(canvas, doc.leftMargin, L1Pos) #Draws line on canvas
        #Grid With Sample Information table categories
        hCats = [["Patient Name","Sample ID","Physician"],
                 ["NIH ID","Sample Date","Study"],
                 ["Sequence ID","Sample Type","Clade"]]
        #Converts categories and data to formatted paragraphs
        data = []
        for i, val in enumerate(hCats):
            data.append([])
            for valu in val:
                hData = metadata.get(valu, "")
                data[i].append(Paragraph('<b>'+valu+':</b> '+str(hData), style=styles['std9']))

        #creates a style for the table
        dataTable = Table(data, colWidths=[252,144,144])
        tblStyle = TableStyle([('FONTSIZE',(0,0),(-1,-1), 9),
                               ('FONT', (0,0), (-1,-1), 'Times'),
                               ('VALIGN',(0,0),(-1,-1),'TOP'),
                               ('BOTTOMPADDING',(0,0),(-1,-1),0),
                               ('TOPPADDING',(0,0),(-1,-1),0),
                               ('LEFTPADDING',(0,0),(-1,-1),0),
                               ('ENDPADDING',(0,0),(-1,-1),0)])
        dataTable.setStyle(tblStyle)
        tW, tH = dataTable.wrap(doc.width, doc.topMargin)
        dataTablePos = L1Pos - tH - 2.5
        dataTable.drawOn(canvas, doc.leftMargin, dataTablePos)

        #Second Horizontal Rule
        hL2 = lineFlow(doc._rightMargin - doc.leftMargin, 0)
        lW2, lH2 = hL1.wrap(doc.width, doc.topMargin)
        l2Pos = dataTablePos - lH2 - 8.2
        hL2.drawOn(canvas, doc.leftMargin, l2Pos)

        #FOOTER
        #Creates data for last line. VISL, Page Number, Report Date
        fTData = [["Virus Isolation and Serology Laboratory"]]
        pStr = "Page {}".format(canvas.getPageNumber())
        fTData[0].append(pStr)
        fTDPara = Paragraph("<b>Report Date: </b>"+ datetime.now().strftime("%m/%d/%y"),
                            style=ParagraphStyle('temp',leading=12,
                                       fontSize=9,fontName='Times',alignment=2))
        fTData[0].append(fTDPara)
        #Uses table to align text Left, Center, Right on one line
        fTable = Table(fTData, colWidths=[doc.width//3,doc.width//3,doc.width//3], rowHeights=[14.6])
        fStyle = TableStyle([('FONTSIZE',(0,0),(-1,-1), 9),
                             ('FONT',(0,0),(-1,-1), 'Times'),
                             ('ALIGN',(0,0),(0,0), 'LEFT'),
                             ('ALIGN',(1,0),(1,0), 'CENTER'),
                             ('ALIGN',(2,0),(2,0), 'RIGHT'),
                             ('LEFTPADDING',(0,0),(-1,-1),0),
                             ('RIGHTPADDING',(0,0),(-1,-1),0)])
        fTable.setStyle(fStyle)
        fTW, fTH = fTable.wrap(doc.width, doc.topMargin)
        fTable.drawOn(canvas, doc.leftMargin, 36)

        #Footer horizontal rule
        fHL = lineFlow(doc._rightMargin - doc.leftMargin, 0)
        flW, flH = fHL.wrap(doc.width,doc.topMargin)
        fHL.drawOn(canvas, doc.leftMargin, 36 + fTH + 5)
        canvas.restoreState()
        return 46 + fTH


    styles = createStylesheet()
    #Creates document template
    print(response)
    doc = StdDocTemplate(response, rightMargin = 36, leftMargin = 36,
                         topMargin = 9.7, bottomMargin = 0, pagesize=letter)
    #Frame for first page has to be smaller than the rest because of footer

    fframe = Frame(doc.leftMargin, 108, doc.width, doc.height - 154.6 - 108,
                  id='normal', leftPadding=0, rightPadding=0, topPadding=0)
    frame = Frame(doc.leftMargin, 60, doc.width, doc.height - 214.6,
                  id='normal', leftPadding=0, rightPadding=0, topPadding=0)
    #Templates for each page
    page = PageTemplate(id='standard', frames = frame, onPage = _header_footer)
    fPage = PageTemplate(id='first', frames = fframe, onPage = _header_footerFP)
    pages=[fPage,page]
    doc.addPageTemplates(pages)
    elements = []

    #Table Style for all four of the drug categories
    resStyle = TableStyle([('FONTSIZE',(0,0),(-1,-1), 8),
                           ('FONT', (0,1), (-1,-1), 'Times'),
                           ('FONT', (0,0), (-1,0), 'TimesBd'),
                           ('GRID', (0,0), (-1,-1), 1, colors.black),
                           ('LEFTPADDING', (0,0), (-1,-1), 2.5),
                           ('TOPPADDING', (0,0), (-1,-1), 2.5),
                           ('BOTTOMPADDING',(0,0),(-1,-1), 0),
                           ('ENDPADDING',(0,0),(-1,-1), 0),
                           ('BOTTOM', (0,0), (-1,-1), 0),
                           ('VALIGN', (0,0), (-1,-1), 'TOP'),
                           ('BACKGROUND', (0,0),(-1,0), colors.burlywood)])
    #Column widths for the drug results tables
    resTblWidths = [86.4,273.6,25.2,28.8,25.2,100.8]
    #List of colors for the "Color" column
    colorArr = [colors.HexColor(0x018523),colors.HexColor(0xffab1a),
                colors.HexColor(0xffab1a),colors.HexColor(0xff6d33),
                colors.HexColor(0xff1100)]
    #Spacer to move tables away from header slightly
    elements.append(ConditionalSpacer(doc.width, 2.2))
    #Used to wrap "mutations" column if it has more than one line of data
    wrapper = TextWrapper(replace_whitespace=True, drop_whitespace=True,width=70,
                          initial_indent='',subsequent_indent='')
    #Titles of each results table
    # json.loads()
    titles = {"NRTI":"Nucleoside Reverse Transcriptase Inhibitors (NRTI)",
              "NNRTI":"Non Nucleoside Reverse Transcriptase Inhibitors (NNRTI)",
              "PR":"Protease Inhibitors (PI)", "IN":"Integrase Inhibitors (INI)"}
    abbrs = ("NRTI", "NNRTI", "PR", "IN")
    for key in abbrs:
        if key in tableData.keys():
            value = tableData[key]
            #Creates table style that uses background colors to populate "Color" column
            subTblStyle = TableStyle(parent=resStyle)
            for row, val in enumerate(value[1:]):
                    subTblStyle.add('BACKGROUND', (4,row+1), (4,row+1), colorArr[val[3]-1])
                    val[1]=wrapper.fill(val[1])
            #Creates title of the table
            elements.append(Paragraph(titles[key], style=styles['tblHead']))
            #Creates actual results table
            elements.append(Table(tableData[key], colWidths=resTblWidths,
                                  style=subTblStyle, repeatRows=1))
            elements.append(ConditionalSpacer(doc.width,6))

    elements.append(Paragraph('<b>Drug resistance algorithm:</b> STANFORD ('+
                              version+')',
                              style=styles['std']))

    elements.append(PageBreak())
    #Spacer used between all elements of the rest of the pdf
    elemSpacer = ConditionalSpacer(doc.width,13)
    #Subtype information header
    elements.append(Paragraph('Subtype Information', style=styles['detHead']))
    elements.append(elemSpacer)
    #Subtype information data
    subtStr = '<b>Subtype and % similarity to closest reference isolate:</b><br/>'
    for value in fileIn['alignedGeneSequences']:
        subtStr+='<b>'+value['gene']['name']+':</b> '+fileIn['bestMatchingSubtype']['displayWithoutDistance']+\
            ' ('+'%.2f'%round(value['matchPcnt'],2)+'%)<br/>'

    elements.append(Paragraph(subtStr,style=styles['detBody']))
    elements.append(elemSpacer)
    #Used to wrap text and create indentation in Interpretation paragraphs
    wrapper = TextWrapper(replace_whitespace=False, drop_whitespace=True,width=168,
                          initial_indent='        ',subsequent_indent='           ')

    for value in fileIn['drugResistance']:
        interpHead = ""
        for category in value['mutationsByTypes']:
            #Creates Line Categories of Mutation list
            if value["gene"]["name"] == "RT" and not category["mutationType"]=="Other":
                interpHead+="<b>"+category["mutationType"]+" Resistance Mutations:</b> "
            else:
                interpHead+="<b>"+value["gene"]["name"]+" "+category["mutationType"]+\
                    " Resistance Mutations:</b> "
            #Adds mutations to list above mutations
            if len(category["mutations"]) == 0:
                interpHead+=' '
            for mutation in category["mutations"]:
                interpHead+=mutation["text"]+", "
            interpHead = interpHead[:-2]
            interpHead+="<br/>"
        #Creates gene header
        elements.append(Paragraph("Drug Resistance: "+value["gene"]["name"],style=styles['detHead']))
        elements.append(elemSpacer)
        #Creates box with list of mutations separated by category
        elements.append(Paragraph(interpHead, style=styles['detBody']))
        elements.append(elemSpacer)
        #Header for interpretation paragraph
        interpCat = "<b>"+value["gene"]["name"]+" Interpretation:</b>\n\n"
        interpTxt = ""
        #Creates unformatted paragraph text
        for category in value['commentsByTypes']:
            interpTxt+="<b>"+category["commentType"]+"</b>\n"
            #Adds comments to paragraph text
            for comment in category["comments"]:
                interpTxt+="&bull; "+comment["text"]+"\n"
            interpTxt+="\n"
        #Uses wrapper to reformat text line by line
        for line in interpTxt.splitlines():
            interpCat+=wrapper.fill(line)+'\n'
        #Use XPreformatted to create paragraph while preserving indentation whitespace
        elements.append(XPreformatted(interpCat+"<br/>", style=styles['detLists']))
        elements.append(KeepTogether(elemSpacer))

    #Disclaimer
    elements.append(Paragraph('VISL Disclaimer', style=styles['detHead']))
    elements.append(ConditionalSpacer(doc.width, 3))
    elements.append(Paragraph('''<b>Each VISL Resistance report is manually reviewed
                                by the VISL lab head before being reported to the
                                NIH. Each report must be used in conjunction with a
                                patient's clinical history, and an advanced grasp of
                                 HIV ARV treatment, in order for it to be useful in
                                 a clinical setting.</b>''', style=styles['detBody']))

    #Creates and saves final document
    doc.build(elements)
    #print(doc.width)
    return response

if __name__ == '__main__':
    main(sys.argv[1],sys.argv[2])
