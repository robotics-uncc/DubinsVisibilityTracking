from typing import Any
import yaml


class ConfigurationFactory:
    __instance = None

    @staticmethod
    def getInstance(file: 'str | None' = None) -> 'dict[str, Any]':
        if ConfigurationFactory.__instance is not None:
            return ConfigurationFactory.__instance
        if file is None:
            return None
        with open(file) as f:
            ConfigurationFactory.__instance = yaml.load(f, yaml.FullLoader)
            return ConfigurationFactory.__instance
