import streamlit as st
from multiapp import MultiApp
from apps import home, mol_vis, reactividad, parametros

app = MultiApp()

st.markdown("" "")

app.add_app("Inicio",home.app)
app.add_app("Visualización molecular",mol_vis.app)
app.add_app("Calculo de parámetros de reactividad",reactividad.app)
app.add_app("Otras propiedades",parametros.app)

app.run()
