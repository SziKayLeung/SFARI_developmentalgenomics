from flask import Flask, render_template, url_for
import pandas as pd
import json
import plotly
import plotly.express as px
import plotly.graph_objects as go


app = Flask(__name__)


posts = [
    {
        'author': 'Complex Disease Epigenetics group',
        'title': 'Dataset',
        'content': 'Long-read data is great!',
        'date_posted': '29 May 2023'

    }

]

@app.route("/")
@app.route("/home")
def home():
    return render_template("home.html", posts=posts)

@app.route("/about")
def about():
    return render_template("about.html", title="About")

@app.route("/stats")
def bar_with_plotly():
    # Students data available in a list of list
    students = [['Trem2', 'Iso1', 'Control', 'F', 40],
                ['Trem2', 'Iso2', 'Control', 'F', 10],
                ['Trem2', 'Iso1', 'AD', 'M', 100],
                ['Trem2', 'Iso1', 'Control', 'M', 10],
                ['Apoe', 'Iso3', 'Control', 'F', 5],
                ['Apoe', 'Iso3', 'AD', 'M', 0]]

    # Convert list to dataframe and assign column values
    df = pd.DataFrame(students,
                      columns=['Gene', 'Isoform', 'Genotype', 'Sex', 'Counts'])

    # Create box plot chart
    fig = px.box(df, x='Isoform', y='Counts', color='Sex', points='all')
    fig.update_layout(title='Box-plot of expression', xaxis_title='Isoform', yaxis_title='Expression', width=800, height=800)

    buttons = [
        dict(
            args=[{'y': [df[df['Gene'] == gene]['Counts']]}],
            label=gene,
            method="update"
        )
        for gene in df['Gene'].unique()
    ]

    fig.update_layout(
        updatemenus=[
            dict(
                buttons=buttons,
                direction="down",
                pad={"r": 10, "t": 10},
                showactive=True,
                x=0.15,  # Adjust the x value to align the updatemenus
                xanchor="left",
                y=1.15,
                yanchor="top"
            ),
        ],
        annotations=[
            dict(
                text="Genes",
                showarrow=False,
                xref="paper",
                yref="paper",
                x=0.05,  # Adjust the x value to align the annotation text
                y=1.12,
                font=dict(size=15)
            )
        ],
        margin=dict(t=200)  # Adjust the top margin value to increase the space between title and plot
    )

    # Convert figure to JSON
    graphJSON = json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)

    # Use render_template to pass graphJSON to HTML
    return render_template('bar.html', graphJSON=graphJSON)


@app.route("/visualisation")
def visualise():
    return render_template("visualisation.html", title="Visualisation")

if __name__ == "__main__":
    app.run(debug=True)
