#!/usr/bin/env python3

from shiny import App, render, ui

app_ui = ui.page_fluid(
    ui.h3("Py-Shiny Jupyter Smoke Test"),
    ui.input_slider("n", "Choose a number", 1, 100, 25),
    ui.output_text_verbatim("out"),
)


def server(input, output, session):
    @output
    @render.text
    def out():
        return f"Selected value: {input.n()}"


app = App(app_ui, server)

