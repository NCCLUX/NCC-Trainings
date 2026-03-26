from datetime import datetime

def on_config(config):
    if config['copyright']:
        config['copyright'] = config['copyright'].format(current_year=datetime.now().year)
    return config
