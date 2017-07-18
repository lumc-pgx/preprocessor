import pandas as pd
from pandasql import sqldf
import collections

from bokeh.plotting import figure, output_file, save
from bokeh.models import HoverTool, Range1d
from bokeh.models.widgets import Div
from bokeh.layouts import Column, Row, ToolbarBox

def parse_data(zmw_data, var_data):
    # join variant and zmw_data
    query = """
        SELECT zmw_data.*, var_data.Pos, var_data.Length, var_data.Type, 
               var_data.RefBP, var_data.AltBP, var_data.QV, var_data.IndelType, var_data.Bases
        FROM zmw_data, var_data
        WHERE zmw_data.ZMW = var_data.ZMW
        """

    if not snakemake.params["show_indel"]:
        query += " and var_data.Type = 'SNP'"

    query += ";"
    
    zmw_with_variants = sqldf(query, locals())
    
    if snakemake.params["show_indel"]:
        zmw_with_variants.loc[zmw_with_variants["IndelType"] == "Deletion", "RefBP"] = \
            zmw_with_variants[zmw_with_variants["IndelType"] == "Deletion"]["Bases"]
        
        zmw_with_variants.loc[zmw_with_variants["IndelType"] == "Insertion", "AltBP"] = \
            zmw_with_variants[zmw_with_variants["IndelType"] == "Insertion"]["Bases"]
    
        zmw_with_variants.loc[zmw_with_variants["Type"] == "INDEL", "Type"] = \
            zmw_with_variants[zmw_with_variants["Type"] == "INDEL"]["IndelType"]
        
    # Get variants with QV above threshold (restrict to only those which hit specified chromosome) 
    high_q = sqldf(
       """
       SELECT *
       FROM zmw_with_variants 
       WHERE QV >= {qv} AND Ref = '{chrom}' and ReadLength >= {length};
       """.format(qv=snakemake.params.min_qv, chrom=snakemake.params.chrom, length=snakemake.params.min_length), locals())
    
    num_molecules = sqldf(
        """
        SELECT COUNT (*) AS mols FROM (SELECT DISTINCT Movie, ZMW FROM high_q);
        """, locals())["mols"][0]
    
    # Count variants
    count_table = sqldf(
        """
        SELECT Pos, Length, Type, RefBP, AltBP, 
               count(*) as Count, sum(NP) as Passes, 
               CAST(count(*) as float) / CAST({} as float) as Freq
        FROM high_q
        GROUP BY Pos, Length, Type, RefBP, AltBP;
        """.format(num_molecules), locals())
    
    return count_table


def make_plot(source_data, title, link_axes=None):
    p = figure(plot_width=800, plot_height=300, logo=None, toolbar_location=None)
    p.title.text = title
    p.xaxis.axis_label = "Position ({})".format(snakemake.params.chrom)
    p.yaxis.axis_label = "Normalized Frequency"
    
    if link_axes is not None:
        p.x_range = link_axes.x_range
        p.y_range = link_axes.y_range
    else:
        p.y_range = Range1d(0, 1.05)
    
    # filter out low frequency
    data = source_data[source_data["Freq"] > snakemake.params.min_freq]
    
    p.x(x="Pos", y="Freq", source=data[data["Type"] == "SNP"],
        size=5, color="navy", alpha=0.75, line_width=1)
    
    p.cross(x="Pos", y="Freq", source=data[data["Type"] == "Insertion"],
            size=7, color="navy", alpha=0.75, line_width=1)
    
    p.circle(x="Pos", y="Freq", source=data[data["Type"] == "Deletion"],
             size=5, line_color="navy", line_alpha=0.75, line_width=1, fill_alpha=0)
    
    hover = HoverTool(tooltips=[
        ("Pos", "@Pos"),
        ("Type", "@Type"),
        ("Ref", "@RefBP"),
        ("Alt", "@AltBP"),
    ])

    p.add_tools(hover)

    return p


def dataframe_from_csv(csv_filenames):
    print(csv_filenames)
    if isinstance(csv_filenames, collections.Iterable) and not isinstance(csv_filenames, str):
        data = pd.read_csv(csv_filenames[0])
        if len(csv_filenames) > 1:
            for filename in csv_filenames[1:]:
                data.append(pd.read_csv(filename))
    else:
        data = pd.read_csv(csv_filenames)
    
    return data.fillna('')


def filter_subreads(data, subreads):
    """
    filter 'data' to only include reads present in 'subreads'
    """
    return sqldf(
        """
        SELECT * FROM data
        WHERE data.ZMW IN (SELECT ZMW FROM subreads) AND
              data.Movie IN (SELECT Movie FROM subreads);
        """, locals())

 
def process(zmw_csv, variant_csv, laa_csv, laa_summary_csv):
    zmw_data = dataframe_from_csv(zmw_csv)
    var_data = dataframe_from_csv(variant_csv)
    
    counts = parse_data(zmw_data, var_data)
    plots = Column(make_plot(counts, "Consensus All"))
    tools = plots.children[0].tools
    
    laa_data = dataframe_from_csv(laa_csv)
    laa_data["Movie"] = laa_data["SubreadId"].map(lambda x: x.split("/")[0])
    laa_data["ZMW"] = laa_data["SubreadId"].map(lambda x: int(x.split("/")[1]))
    
    laa_summary = dataframe_from_csv(laa_summary_csv)
    laa_summary = laa_summary[(laa_summary["NoiseSequence"] == 0) & (laa_summary["IsDuplicate"] == 0) & (laa_summary["IsChimera"] == 0)]
    good_amplicons = list(laa_summary["FastaName"])
    
    for amplicon in good_amplicons:
        # pick the subreads for this amplicon
        subreads = laa_data[laa_data[amplicon] == 1]
        
        # select the corresponding subreads from the zmw and variant data
        zmws = filter_subreads(zmw_data, subreads)
        variants = filter_subreads(var_data, subreads)
        
        # make a plot
        counts = parse_data(zmws, variants)
        plt = make_plot(counts, amplicon, plots.children[0])
        tools += plt.tools
        plots.children.append(plt)
    
    subreads_remaining = laa_data[laa_data[good_amplicons].sum(axis=1) == 0]
    # select the remaining subreads from the zmw and variant data
    zmws = filter_subreads(zmw_data, subreads_remaining)
    variants = filter_subreads(var_data, subreads_remaining)
    
    # make a plot
    counts = parse_data(zmws, variants)
    plt = make_plot(counts, "Consensus Remaining", plots.children[0])
    tools += plt.tools
    plots.children.append(plt)
    
    # merge the toolbars
    toolbar = ToolbarBox(
        logo=None,
        merge_tools=True,
        toolbar_location="right",
        tools=tools,
    )
    
    return Row(plots, toolbar)


output_file(snakemake.output[0], title=snakemake.wildcards.barcode)

try:
    plots = process(
        zmw_csv = snakemake.input.ccs_zmws,
        variant_csv = snakemake.input.ccs_variants,
        laa_csv = snakemake.input.laa_subreads,
        laa_summary_csv = snakemake.input.laa_summary
    )
    
    save(plots)
except pd.errors.EmptyDataError:
    save(Div(text="No Data"))
