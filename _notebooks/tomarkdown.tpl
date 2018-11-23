{% extends 'markdown.tpl' %}


{% block stream %}
{{"[OUTPUT]" | indent}} 

{{ output.text | indent}}
{% endblock stream %}

{% block data_text scoped %}
{{"[OUTPUT]" | indent}} 

{{ output.data['text/plain'] | indent }}
{% endblock data_text %}



