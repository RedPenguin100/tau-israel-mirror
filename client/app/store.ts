import { action, atom, map } from "nanostores";

export const isLoading = atom(false);

export function setIsLoading(newLoadingState: boolean) {
  isLoading.set(newLoadingState);
}

export type Annotation = {
  name: string;
  start: number;
  end: number;
  direction: 1 | -1;
};

export type PlasmidMap = {
  name: string;
  sequence: string;
  annotations: Annotation[];
};

export const plasmidMap = map<PlasmidMap>({
  name: "",
  sequence: "",
  annotations: [],
});

export const setPlasmidMap = action(
  plasmidMap,
  "setPlasmidMap",
  (store, { name, sequence, annotations }: PlasmidMap) => {
    store.setKey("name", name);
    store.setKey("sequence", sequence);
    store.setKey("annotations", annotations);
  }
);

export const resetPlasmidMap = action(
  plasmidMap,
  "resetPlasmidMap",
  (store) => {
    store.setKey("name", "");
    store.setKey("sequence", "");
    store.setKey("annotations", []);
  }
);
